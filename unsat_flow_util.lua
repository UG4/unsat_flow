-- author: Niklas Conen

-- utility functions for unsaturated density flow equations
-- problems are defined using config files

ug_load_script("ug_util.lua")

-- namespace creation
util.unsat = util.unsat or {}

local json = util.json


-- object oriented clustering of the Approximation Space, Element Discretisation and
-- Domain Discretisation
ProblemDisc = {}

function ProblemDisc:new(problemDesc, dom)
    assert(problemDesc, "No Input defined")
    problemObj = {}
    setmetatable(problemObj, self)
    self.__index = self

    self.problem = problemDesc
    self.domain = dom
    self.cmp = problemDesc.flow.cmp

    self.gravity = ConstUserVector(0.0)
    self.gravity:set_entry(problemDesc.domain.dim-1, problemDesc.flow.gravity)

    self.modelMap = nil
    if problemDesc.parameter ~= nil then
        print (problemDesc.parameter)
        self.modelMap = ProblemDisc:CreateModelMap(problemDesc.parameter)
    end

    return problemObj
end

function ProblemDisc:CreateApproxSpace()
    -- documentation for ug function ApproximationSpace: line 277 in /ug4/ugcore/ugbase/lib_disc/function_spaces/approximation_space.h
    self.approxSpace = ApproximationSpace(self.domain)
    self.approxSpace:add_fct(self.cmp[1], "Lagrange", 1)
    self.approxSpace:add_fct(self.cmp[2], "Lagrange", 1)
    self.approxSpace:init_levels()
    self.approxSpace:init_top_surface()
    self.approxSpace:print_statistic()
    return self.approxSpace
end

-- Element discretisation
function ProblemDisc:CreateElemDisc(subdom, medium)

  local myUpwind = FullUpwind()
    -- Creates the elememt discretisation for a given medium
	local elemDisc = {}
   
    -- Eqns. for flow and transport.
	elemDisc["flow"] = ConvectionDiffusion(self.cmp[1], subdom, "fv1")      
    elemDisc["transport"] = ConvectionDiffusion(self.cmp[2], subdom, "fv1") 

    -- Capillary pressure (auxiliary).
    local capillary = ScaleAddLinkerNumber()
    capillary:add(-1.0, elemDisc["flow"]:value())

    -- Local viscosity = self.mu
    local viscosity = self:viscosity()
    if self.problem.flow.viscosity.type ~= "const" then
       viscosity:set_input(0, elemDisc["transport"]:value())
    end


    local density = self:density()
    density:set_input(0, elemDisc["transport"]:value())

    local porosity = medium.porosity
    local permeability = medium.permeability

    -- the vanGenuchten model is calculated using the Richards Plugin
    local conductivity = ProblemDisc:conductivity(medium.conductivity.value, capillary)
    local saturation = ProblemDisc:saturation(medium.saturation.value, capillary)

    -- hydraulic conductivity in saturated medium
    -- K_s = k / mu
    -- for unsaturated flow: 
    -- multiplied by relative hydraulic conductivity k(p)
    local khyd = ScaleAddLinkerMatrix()
    khyd:add(conductivity, permeability)

    -- Darcy Velocity
    -- $\vec q := -k*k(p)/mu (\grad p - \rho \vec g)$
    local DarcyVelocity = DarcyVelocityLinker()
    DarcyVelocity:set_permeability(khyd)  
    DarcyVelocity:set_viscosity(viscosity)
    DarcyVelocity:set_pressure_gradient(elemDisc["flow"]:gradient())
    DarcyVelocity:set_density(density)
    DarcyVelocity:set_gravity(self.gravity)
	-- print("Darcy Velocity created.")

    -- 
   

    -- molecular diffusion
    local transportDiffTensor = ScaleAddLinkerMatrix()
    transportDiffTensor:add(porosity*density, medium.diffusion)

    if (medium.alphaT) then 
        local dispersion = BearScheidegger(DarcyVelocity, medium.alphaT, medium.alphaL)
        transportDiffTensor:add(density, dispersion)
    end

	-----------------------------------------
	-- Equation [1]
	-----------------------------------------
	-- flow equation
	-- $\partial_t (\Phi \rho_w S_w) 
	--		+ \nabla \cdot [\rho_w \vec{v}_w] = \rho_w \Gamma_w$

	-- fluid storage: \Phi \rho_w S_w
    local fluidStorage = ScaleAddLinkerNumber()
	-- flux \rho_w \vec{v}_w
    local fluidFlux = ScaleAddLinkerVector()
    
    fluidStorage:add(porosity*density, saturation) -- INSERT: saturation  or 1.0
    fluidFlux:add(density, DarcyVelocity)

    if self.problem.flow.boussinesq then
        -- boussinesq fluid storage: \Phi \rho_w0 S_w
        local storageB = ScaleAddLinkerNumber()
        -- boussinesq flux \rho_w0 \vec{v}_w
        local fluidFluxB = ScaleAddLinkerVector()

        storageB:add(porosity*self.problem.flow.density.min, saturation)
        fluidFluxB:add(self.problem.flow.density.min, DarcyVelocity)

        elemDisc["flow"]:set_mass(storageB)
        elemDisc["flow"]:set_flux(fluidFluxB)

    else
       

        elemDisc["flow"]:set_mass(fluidStorage)
        elemDisc["flow"]:set_flux(fluidFlux)
    end


	elemDisc["flow"]:set_mass_scale(0.0)
    elemDisc["flow"]:set_diffusion(0.0)

	-----------------------------------------
	-- Equation [2]
	-----------------------------------------
	-- transport convection and diffusion => transport equation
	-- $\partial_t (\Phi \rho_w S_w \omega) 
	-- 		+ \nabla \cdot [\rho_w \omega \vec{v}_w - \rho_w D \nabla \omega] 
	--			= \rho_w \omega \Gamma_w$

    elemDisc["transport"]:set_mass(0.0)
    elemDisc["transport"]:set_mass_scale(fluidStorage) -- * w
    elemDisc["transport"]:set_velocity(fluidFlux)
    elemDisc["transport"]:set_diffusion(transportDiffTensor)
    elemDisc["transport"]:set_upwind(myUpwind)


   
 --   ScaleAddLinkerNumber()
  --  capillary:add(-1.0, elemDisc["flow"]:value())
    
   --[[ if medium.conductivity.type == "exp" then
        conductivity:set_input(0, capillary)

    elseif medium.conductivity.type == "vanGenuchten" then
        conductivity:set_capillary(capillary)
    end

    if medium.saturation.type == "exp" then
        saturation:set_input(0, capillary)

    elseif medium.saturation.type == "vanGenuchten" then
        saturation:set_capillary(capillary)
    end
    --]]

    print("Created Element Discretisation for Subset ", subdom)


    -- Preparations for IO.
    local si = self.domain:subset_handler():get_subset_index(subdom)
    self.CompositeCapillary:add(si, capillary)
    self.CompositePermeability:add(si, conductivity)
    self.CompositeSaturation:add(si, saturation)
    self.CompositeDarcyVelocity:add(si, DarcyVelocity)
    self.CompositeFlux:add(si, fluidFlux)

    return elemDisc
end


function ProblemDisc:CreateDomainDisc(approxSpace)
    self.domainDisc = DomainDiscretization(approxSpace)

    local conc = GridFunctionNumberData(self.u, "c")

    -- density grid function
    -- self.rho = self:density()
    --self.rho:set_input(0, conc)

    self.mu = self:viscosity()
    if self.problem.flow.viscosity.type ~= "const" then
       -- self.mu:set_input(0, conc)
    end

    -- self.grad_p = GridFunctionGradientData(self.u, "p")

    --self.capillary = ScaleAddLinkerNumber()
    --self.capillary:add(-1.0, GridFunctionNumberData(self.u, "p"))

    if self.problem.boussinesq then
        print("using boussineq approximation")
    end
    
    self.CompositeCapillary = CompositeUserNumber(true)
    self.CompositePermeability = CompositeUserNumber(false)
    self.CompositeSaturation = CompositeUserNumber(false)
    self.CompositeDarcyVelocity = CompositeUserVector(false)
    self.CompositeFlux = CompositeUserVector(false)

    for i,medium in ipairs(self.problem.medium) do -- for all media
        local elemDisc = nil

        for j, subset in ipairs(medium.subsets) do
            elemDisc = self:CreateElemDisc(subset, medium)

            self.domainDisc:add(elemDisc["flow"])
            self.domainDisc:add(elemDisc["transport"])
        end
    end

    local dim = self.problem.domain.dim
    -- Sources and sinks
    if (self.problem.sources) then
    for i,mySource in ipairs(self.problem.sources) do -- for all sources
        
        print(mySource)
       
        -- Read Vector.
        local mycoord = Vec()
        mycoord:set_coord(0, mySource.coord[1])
        mycoord:set_coord(1, mySource.coord[2])
        if (dim==3) then mycoord:set_coord(2, mySource.coord[3]) end 
        -- Add primary.
        local elemDisc = DiracSourceDisc(mySource.cmp, mySource.subset)
        elemDisc:add_source(mySource.strength, mycoord)
        self.domainDisc:add(elemDisc)
        
        print("Added point source for " .. mySource.cmp.. " at ");
        print(mySource.strength)
        print(mySource.coord)

    
        -- Add source for substances
        print (mySource.substances)
        for j,mySub in ipairs(mySource.substances) do -- for all components
            local elemDiscSub = DiracSourceDisc(mySub.cmp, mySource.subset)
            elemDiscSub:add_transport_sink(mySource.strength)
            self.domainDisc:add(elemDiscSub)
        end

       
    end
    end

    -- Sinks
    if (self.problem.sinks) then 
    for index,mySink in ipairs(self.problem.sinks) do -- for all sources
    
        print (mySink)
        print(mySink.strength)

        -- Create vector.
        local mycoord = Vec()
        mycoord:set_coord(0, mySink.coord[1])
        mycoord:set_coord(1, mySink.coord[2])
        if (dim==3) then mycoord:set_coord(2, mySource.coord[3]) end 

        -- Add primary.
        local elemDisc = DiracSourceDisc(mySink.cmp, mySink.subset)
        elemDisc:add_transport_sink(mySink.strength)
        self.domainDisc:add(elemDisc)
        print("Added sink for" .. mySource.cmp)

        -- Add source for substances

      end
    end
    
    -- Create Boundary Conditions
    local dirichletBnd = nil
    local neumannBnd = {} 
    for i, v in ipairs(self.problem.boundary) do
       
        if v.type == "dirichlet" then
            -- Dirichtlet boundary 
            dirichletBnd = dirichletBnd or DirichletBoundary()
            dirichletBnd:add(v.value, v.cmp, v.bnd)
            print("Added Dirichlet Boundary with value " .. v.value .. " for " .. v.cmp .. " on subset " .. v.bnd)
        end
        
        if v.type == "flux" then
            -- Neumann-type
            neumannBnd[v.cmp] =  neumannBnd[v.cmp] or NeumannBoundary(v.cmp, "fv1")
            neumannBnd[v.cmp]:add(v.value, v.bnd, v.inner)
            print("Added Neumann Boundary with value " .. v.value .. " for " .. v.cmp .. " on subset " .. v.bnd)
        end
 
    end
    
    if (dirichletBnd) then  self.domainDisc:add(dirichletBnd) end
    if (neumannBnd["p"]) then  self.domainDisc:add(neumannBnd["p"]) end
    if (neumannBnd["c"]) then  self.domainDisc:add(neumannBnd["c"]) end


    local constraints = OneSideP1Constraints()
    self.domainDisc:add(constraints)
    
    print("Created Domain Discretisation")
    return self.domainDisc
end


function ProblemDisc:CreateVTKOutput()
    -- generating the vtk output
    -- for future reference; util.Balance()
    for i, v in ipairs(self.problem.output.data) do
        -- concentration and pressure
        if v == "p" or v == "c" then
           self.vtk:select_nodal(GridFunctionNumberData(self.u, v), v)
        -- density
        elseif v == "rho" then
            --self.vtk:select_element(self.rho, v)
        -- viscosity
        elseif v == "mu" and self.problem.flow.viscosity.type ~= "const" then
            --self.vtk:select_element(self.mu, v)
        -- conductivity
        elseif v == "k" then
            self.vtk:select_element(self.CompositePermeability, v)
        -- saturation
        elseif v == "s" then
            self.vtk:select_element(self.CompositeSaturation, v)
        -- darcy velocity
        elseif v == "q" then
            self.vtk:select_element(self.CompositeDarcyVelocity, v)
        -- transport equation flux
        elseif v == "f" then
            self.vtk:select_element(self.CompositeFlux, v)
        end
        
    end
    self.vtk:select(self.CompositeCapillary, "pc")
end


function ProblemDisc:CreateModelMap(paramDesc)
    local modelMap = {}
    local mfactory = RichardsModelFactory()
          
    for i, medium in ipairs(paramDesc) do
        if medium.type == "vanGenuchten" then
            modelMap[medium.uid] = mfactory:create_van_genuchten(json.encode(medium))
            print(medium.uid..":v->"..modelMap[medium.uid]:config_string())
        elseif medium.type == "exp" then
            modelMap[medium.uid] = mfactory:create_exponential(json.encode(medium))
            print(medium.uid..":e->"..modelMap[medium.uid]:config_string())
        elseif medium.type == "const" then
            modelMap[medium.uid] = medium.value
            print(medium.uid..":c->"..medium.value)
        end
    end
    return modelMap
end


function ProblemDisc:SetInitialData()
    for i, initial in ipairs(self.problem.initial) do
        print("Added Initial Value " .. initial.cmp .. " = " .. initial.value)
        Interpolate(initial.value, self.u, initial.cmp)
    end
end


function ProblemDisc:density(densDesc)   
    local p_w = self.problem.flow.density.min
    local p_s = self.problem.flow.density.max
    local w_max = self.problem.flow.density.w_max
    local density = nil
    
    -- p_s: max density -> saline
    -- p_w: min density -> water
    -- w: mass fraction
    if self.problem.flow.density.type == "linear" then
        -- linear density function
        function DensityFct(w)
            return p_w + (p_s - p_w) * (w/w_max)
        end

        function dwDensityFct(w)
            return (p_s - p_w) / w_max
        end

        density = LuaUserFunctionNumber("DensityFct", 1);
        density:set_deriv(0, "dwDensityFct");

    elseif self.problem.flow.density.type == "exp" then
        -- exponential density function
        function DensityFct(w)
            return p_w * (p_s / p_w)^(w/w_max)
        end

        function dwDensityFct(w)
            return (p_w * (p_s / p_w)^(w/w_max) * math.log(p_s / p_w)) / w_max
        end

        density = LuaUserFunctionNumber("DensityFct", 1);
        density:set_deriv(0, "dwDensityFct");

    elseif self.problem.flow.density.type == "ideal" then
        -- ideal density function
        function DensityFct(w)
            return 1/(1/p_w + (1/p_s - 1/p_w) * (w/w_max))
        end

        function dwDensityFct(w)
            return (w_max * p_s * (w_max * p_s - w^2))/(w_max * p_s + w * (-1.0 * p_s + w))^2
        end

        density = LuaUserFunctionNumber("DensityFct", 1);
        density:set_deriv(0, "dwDensityFct");
    end

    return density
end

function ProblemDisc:conductivity(condID, pcapillary)
    local conductivity = nil
    local values = self.modelMap[condID]

    local factory = RichardsUserDataFactory(pcapillary)

    if type(values) == "userdata" then
        -- works for exponential and van Genuchten.
        print(values)
        conductivity = factory:create_conductivity(values)
        print(conductivity)
    elseif type(values) == "table" then
        -- TODO: Deprecated. 
        if values.type == "exp" then


            
            local alpha = values.alpha
            local thetaS = values.thetaS
            local thetaR = values.thetaR
            function ConductivityFct(p)
                if p < 0 then
                    return 1
                else
                    local temp = thetaR + (thetaS - thetaR) * math.exp(alpha * p)
                    if temp > 1 then
                        return 1
                    elseif temp < 0 then
                        return 0
                    else 
                        return temp
                    end
                end
            end

            function dpConductivityFct(p)
                if p < 0 then
                    return 0
                else
                    return (thetaS - thetaR) * alpha * math.exp(alpha * p)
                end
            end

            conductivity = LuaUserFunctionNumber("ConductivityFct", 1)
            conductivity:set_deriv(0, "dpConductivityFct")
            conductivity:set_input(0, pcapillary)
        end
    end
    return conductivity
end


function ProblemDisc:saturation(satID, pcapillary)
    local saturation = nil
    local values = self.modelMap[satID]
    local factory = RichardsUserDataFactory(pcapillary)

    if type(values) == "userdata" then -- works for exponential and van Genuchten.
        print(values)
        saturation = factory:create_saturation(values)
        print(saturation)
    elseif type(values) == "table" then
        if values.type == "exp" then
            local alpha = values.alpha
            local K_sat = values.K_sat

            function SaturationFct(p)
                if p < 0 then
                    return 1
                else
                    return K_sat * math.exp(alpha * p)
                end
            end

            function dpSaturationFct(p)
                if p < 0 then
                    return 0
                else
                    return K_sat * alpha * math.exp(alpha * p)
                end
            end

            saturation = LuaUserFunctionNumber("SaturationFct", 1)
            saturation:set_deriv(0, "dpSaturationFct")
            conductivity:set_input(0, pcapillary)
        end   
    end

    return saturation
end


function ProblemDisc:viscosity()
    local viscosity = nil
    local mu0 = nil
    if self.problem.flow.viscosity.type == "const" then
        viscosity = self.problem.flow.viscosity.mu0

    elseif self.problem.flow.viscosity.type == "real" then
        mu0 = self.problem.flow.viscosity.mu0

        function ViscosityFct(w)
            return mu0*(1 + 1.85*w - 4.1*w^2 + 44.5 * w^3)
        end

        function dwViscosityFct(w)
            return mu0*(1.85 - 8.2*w + 133.5*w^2)
        end

        viscosity = LuaUserFunctionNumber("ViscosityFct", 1);
        viscosity:set_deriv(0, "dwViscosityFct");
    end

    return viscosity
end
