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
    self.approxSpace:add_fct("p", "Lagrange", 1)
    self.approxSpace:add_fct("c", "Lagrange", 1)
    self.approxSpace:init_levels()
    self.approxSpace:init_top_surface()
    self.approxSpace:print_statistic()
    return self.approxSpace
end

-- Element discretisation
function ProblemDisc:CreateElemDisc(subdom, medium)
    -- Creates the elememt discretisation for a given medium
	local elemDisc = {}
<<<<<<< HEAD
    -- flow equation
	elemDisc["flow"] = ConvectionDiffusion("p", subdom, "fv1")
	-- transport equation
    elemDisc["transport"] = ConvectionDiffusion("c", subdom, "fv1")

=======
   
    -- Eqns. for flow and transport.
	elemDisc["flow"] = ConvectionDiffusion(self.cmp[1], subdom, "fv1")      
    elemDisc["transport"] = ConvectionDiffusion(self.cmp[2], subdom, "fv1") 

    -- Capillary pressure (auxiliary).
    local capillary = ScaleAddLinkerNumber()
    capillary:add(-1.0, elemDisc["flow"]:value())

    -- Local viscosity = self.mu
>>>>>>> feature-arne
    local viscosity = self:viscosity()
    if self.problem.flow.viscosity.type ~= "const" then
        viscosity:set_input(0, elemDisc["transport"]:value())
    end
    self.mu = viscosity

    local density = self:density()
<<<<<<< HEAD
    if self.problem.flow.density.type ~= "const" then
        density:set_input(0, elemDisc["transport"]:value())
    end
    self.rho = density

    local porosity = nil    -- phi
    if type(medium.porosity) == "string" then
        for i, param in ipairs(self.problem.parameter) do
            if param.uid == medium.porosity then
                porosity = param.thetaS     -- the porosity is equal to the saturated water content
            end
        end
    else
        porosity = medium.porosity
    end

    -- the vanGenuchten model is calculated using the Richards Plugin
    -- capillary pressure: air pressure is set to 0
    -- => p_c = - p_w
    local capillary = -1.0*elemDisc["flow"]:value()
    local RichardsUserData = RichardsUserDataFactory(capillary)
    conductivity = RichardsUserData:create_conductivity(self.modelMap[medium.conductivity.value])
    saturation = RichardsUserData:create_saturation(self.modelMap[medium.saturation.value])
=======
    density:set_input(0, elemDisc["transport"]:value())

    local porosity = medium.porosity
    local permeability = medium.permeability

    -- the vanGenuchten model is calculated using the Richards Plugin
    local conductivity = ProblemDisc:conductivity(medium.conductivity.value, capillary)
    local saturation = ProblemDisc:saturation(medium.saturation.value, capillary)
>>>>>>> feature-arne

    -- hydraulic conductivity in saturated medium is given by darcys law
    -- k_f = K*rho*g / mu
    -- for unsaturated flow:
    -- multiplied by relative hydraulic conductivity K(S)
    -- k_f,r(S) * k_f =  K * rho * g / mu
    Ksat = nil
    for i, param in ipairs(self.problem.parameter) do
        if param.uid == medium.conductivity.value then
            Ksat = param.Ksat
        end
    end

    local r_scale = -86400 -- scaling for time unit in K_sat
    local permeability = ScaleAddLinkerMatrix()
    -- needs scaling depending on hydraulic conductivities scale in the Richardsplugin
    permeability:add(conductivity, self.problem.flow.viscosity.mu0/(r_scale*self.problem.flow.density.min*self.problem.flow.gravity))
    print("saturated permeability: "..(Ksat*self.problem.flow.viscosity.mu0/(r_scale*self.problem.flow.density.min*self.problem.flow.gravity)))

    -- Darcy Velocity
    -- $\vec q := -k*k(p)/mu (\grad p - \rho \vec g)$
    local DarcyVelocity = DarcyVelocityLinker()
    --DarcyVelocity:set_viscosity(self.problem.flow.density.min*self.problem.flow.gravity*r_scale)
    DarcyVelocity:set_viscosity(viscosity)
    DarcyVelocity:set_permeability(permeability)
    DarcyVelocity:set_pressure_gradient(elemDisc["flow"]:gradient())
    DarcyVelocity:set_density(density)
    DarcyVelocity:set_gravity(self.gravity)

<<<<<<< HEAD
    local volufrac = ScaleAddLinkerNumber()
    volufrac:add(porosity, saturation)
=======
    -- 
   

    -- molecular diffusion
    local transportDiffTensor = ScaleAddLinkerMatrix()
    transportDiffTensor:add(porosity*density, medium.diffusion)

    if (medium.alphaT) then 
        local dispersion = BearScheidegger(DarcyVelocity, medium.alphaT, medium.alphaL)
        transportDiffTensor:add(density, dispersion)
    end
>>>>>>> feature-arne

	-----------------------------------------
	-- Equation [1]
	-----------------------------------------
	-- flow equation
	-- $\partial_t (\Phi S_w \rho_w)
	--		+ \nabla \cdot (\rho_w \vec{v}_w) = \rho_w \Gamma_w$

<<<<<<< HEAD
	-- fluid storage: \Phi S_w \rho_w
    local storage = ScaleAddLinkerNumber()
    storage:add(volufrac, density)
	-- flux of the fluid phase \rho_w \vec{v}_w
    local fluidFlux = ScaleAddLinkerVector()
=======
	-- fluid storage: \Phi \rho_w S_w
    local fluidStorage = ScaleAddLinkerNumber()
	-- flux \rho_w \vec{v}_w
    local fluidFlux = ScaleAddLinkerVector()
    
    fluidStorage:add(porosity*density, saturation) -- INSERT: saturation  or 1.0
>>>>>>> feature-arne
    fluidFlux:add(density, DarcyVelocity)

    -- oberbeck-boussinesq approximation
    if self.problem.flow.boussinesq then
        -- boussinesq fluid storage: \Phi \rho_w0 S_w
        local storageB = ScaleAddLinkerNumber()
        -- boussinesq flux \rho_w0 \vec{v}_w
        local fluidFluxB = ScaleAddLinkerVector()

        storageB:add(self.problem.flow.density.min, volufrac)
        fluidFluxB:add(self.problem.flow.density.min, DarcyVelocity)

        elemDisc["flow"]:set_mass(storageB)
        elemDisc["flow"]:set_flux(fluidFluxB)

    else
<<<<<<< HEAD
        elemDisc["flow"]:set_mass(storage)
        elemDisc["flow"]:set_flux(fluidFlux)
    end
    elemDisc["flow"]:set_mass_scale(0.0)
=======
       

        elemDisc["flow"]:set_mass(fluidStorage)
        elemDisc["flow"]:set_flux(fluidFlux)
    end


	elemDisc["flow"]:set_mass_scale(0.0)
    elemDisc["flow"]:set_diffusion(0.0)
>>>>>>> feature-arne

	-----------------------------------------
	-- Equation [2]
	-----------------------------------------
	-- transport convection and diffusion => transport equation
	-- $\partial_t (\Phi \rho_w S_w \omega)
	-- 		+ \nabla \cdot (\rho_w \omega q - \Phi S_w \rho_w,min D \nabla \omega)
	--			= \phi S_w \rho_w \Gamma_w$

<<<<<<< HEAD
    -- diffusive flux: \Phi S_w \rho_w D \nabla \omega
    local diffusion = ScaleAddLinkerMatrix()
    diffusion:add(volufrac, self.problem.flow.diffusion*self.problem.flow.density.min)

    -- advective flux is equal to fluidFlux * mass fraction
    elemDisc["transport"]:set_mass_scale(storage)
    elemDisc["transport"]:set_velocity(fluidFlux)
    elemDisc["transport"]:set_diffusion(diffusion)
    -- set upwind scheme to FullUpwind, if no upwind scheme is defined
    if self.problem.flow.upwind == "full" or self.problem.flow.upwind == nil then
        elemDisc["transport"]:set_upwind(FullUpwind())
    elseif self.problem.flow.upwind == "partial" then
        elemDisc["transport"]:set_upwind(PartialUpwind())
=======
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
>>>>>>> feature-arne
    end
    --]]

    print("Created Element Discretisation for Subset ", subdom)

    -- Preparations for IO.
    local advFlux = ScaleAddLinkerVector()
    c_value = GridFunctionNumberData(self.u, "c")
    advFlux:add(c_value, fluidFlux)

    local gradC = GridFunctionGradientData(self.u, "c")
    local difFlux = ScaleAddLinkerVector()
    difFlux:add(volufrac, self.problem.flow.density.min*self.problem.flow.diffusion*gradC)

    local si = self.domain:subset_handler():get_subset_index(subdom)
    self.CompositeCapillary:add(si, capillary)
    self.CompositeConductivity:add(si, conductivity)
    self.CompositeSaturation:add(si, saturation)
    self.CompositeDarcyVelocity:add(si, DarcyVelocity)
    self.CompositeFluidFlux:add(si, fluidFlux)
    self.CompositeTransportFlux:add(si, (advFlux - difFlux))
    self.CompositeAdvectiveFlux:add(si, advFlux)
    self.CompositeDiffusiveFlux:add(si, difFlux)
    self.CompositeVoluFrac:add(si, volufrac)
    self.CompositeSaltMass:add(si, volufrac*density*elemDisc["transport"]:value())

    return elemDisc
end


function ProblemDisc:CreateDomainDisc(approxSpace)
    self.domainDisc = DomainDiscretization(approxSpace)

    if self.problem.boussinesq then
        print("using boussineq approximation")
    end

    self.CompositeCapillary = CompositeUserNumber(true)
    self.CompositeConductivity = CompositeUserNumber(false)
    self.CompositeSaturation = CompositeUserNumber(false)
    self.CompositeDarcyVelocity = CompositeUserVector(false)
    self.CompositeFluidFlux = CompositeUserVector(false)
    self.CompositeTransportFlux = CompositeUserVector(false)
    self.CompositeAdvectiveFlux = CompositeUserVector(false)
    self.CompositeDiffusiveFlux = CompositeUserVector(false)
    self.CompositeVoluFrac = CompositeUserNumber(false)
    self.CompositeSaltMass = CompositeUserNumber(false)

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
            neumannBnd[v.cmp] = neumannBnd[v.cmp] or NeumannBoundary(v.cmp, "fv1")
            neumannBnd[v.cmp]:add(v.value, v.bnd, v.inner)
            print("Added Neumann Boundary with value " .. v.value .. " for " .. v.cmp .. " on subset " .. v.bnd)
        end

    end

    if (dirichletBnd) then  self.domainDisc:add(dirichletBnd) end
    if (neumannBnd["p"]) then  self.domainDisc:add(neumannBnd["p"]) end
    if (neumannBnd["c"]) then  self.domainDisc:add(neumannBnd["c"]) end

<<<<<<< HEAD
    -- Adding Dirac Sources
    if self.problem.sources ~= nil then
        for i, v in ipairs(self.problem.sources) do

            local source = DiracSourceDisc(v.cmp, v.subset)

            local location = nil
            if self.problem.domain.dim == 2 and v.x and v.y then
                location = Vec2d(v.x, v.y)
            elseif self.problem.domain.dim == 3 and v.x and v.y and v.z then
                location = Vec3d(v.x, v.y, v.z)
            else
                print("source coordinates not given")
            end

            if v.value and type(v.value) ~= "string" then source:add_source(v.value, location)
            elseif v.value then source:add_source(v.value)
            elseif v.transport then source:add_transport_sink(v.transport)
            end
            self.domainDisc:add(source)

            if v.value ~= nil then
                print("Added DiracPointSource with value " .. v.value .. " on subset " .. v.subset)
            else
                print("Added DiracPointSource with transport value of" .. v.transport .. "on subset " .. v.subset)
            end
        end
    end

=======

    local constraints = OneSideP1Constraints()
    self.domainDisc:add(constraints)
    
>>>>>>> feature-arne
    print("Created Domain Discretisation")
    return self.domainDisc
end


function ProblemDisc:CreateVTKOutput()
    -- generating the vtk output
    -- for future reference; util.Balance()
    for i, v in ipairs(self.problem.output.data) do
        -- concentration and pressure
        if v == "p" or v == "all" then
           self.vtk:select_nodal(GridFunctionNumberData(self.u, v), v)
        -- concentration gradient
        elseif v == "gradp" or v == "all" then
            self.vtk:select_element(GridFunctionGradientData(self.u, "p"), v)
        -- concentration
        elseif v == "c" or v == "all" then
            self.vtk:select_nodal(GridFunctionNumberData(self.u, v), v)
        -- concentration gradient
        elseif v == "gradc" or v == "all" then
            self.vtk:select_element(GridFunctionGradientData(self.u, "c"), v)
        -- density
        elseif (v == "rho" or v == "all") and self.problem.flow.density.type ~= "const" then
            self.vtk:select_element(self.rho, v)
        -- viscosity
        elseif (v == "mu" or v == "all") and self.problem.flow.viscosity.type ~= "const" then
            self.vtk:select_element(self.mu, v)
        -- relative conductivity
        elseif v == "kr" or v == "all" then
            self.vtk:select_element(self.CompositeConductivity, v)
        -- saturation
        elseif v == "s" or v == "all" then
            self.vtk:select_element(self.CompositeSaturation, v)
        -- darcy velocity
        elseif v == "q" or v == "all" then
            self.vtk:select_element(self.CompositeDarcyVelocity, v)
        -- flow equation flux
        elseif v == "ff" or v == "all" then
            self.vtk:select_element(self.CompositeFluidFlux, v)
        -- transport equation flux
        elseif v == "tf" or v == "all" then
            self.vtk:select_element(self.CompositeTransportFlux, v)
        -- advective flux
        elseif v == "af" or v == "all" then
            self.vtk:select_element(self.CompositeAdvectiveFlux, v)
        -- diffusive flux
        elseif v == "df" or v == "all" then
            self.vtk:select_element(self.CompositeDiffusiveFlux, v)
        -- capillary pressure
        elseif v == "pc" or v == "all" then
            self.vtk:select(self.CompositeCapillary, v)
        -- volumetric fraction
        elseif v == "vf" or v == "all" then
            self.vtk:select(self.CompositeVoluFrac, v)
        end
    end

end


function ProblemDisc:CreateModelMap(paramDesc)
    local modelMap = {}
    local mfactory = RichardsModelFactory()
          
    for i, medium in ipairs(paramDesc) do
        if medium.type == "vanGenuchten" then
<<<<<<< HEAD
            medium.alpha= medium.alpha / (-1.0 * self.problem.flow.gravity * self.problem.flow.density.min)
            modelMap[medium.uid] = RichardsModelFactory():create_van_genuchten(json.encode(medium))
        elseif medium.type == "const" then
            modelMap[medium.uid] = medium.value
=======
            modelMap[medium.uid] = mfactory:create_van_genuchten(json.encode(medium))
            print(medium.uid..":v->"..modelMap[medium.uid]:config_string())
        elseif medium.type == "exp" then
            modelMap[medium.uid] = mfactory:create_exponential(json.encode(medium))
            print(medium.uid..":e->"..modelMap[medium.uid]:config_string())
        elseif medium.type == "const" then
            modelMap[medium.uid] = medium.value
            print(medium.uid..":c->"..medium.value)
>>>>>>> feature-arne
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
    local density = nil

    -- p_s: max density -> saline
    -- p_w: min density -> water
    -- w: mass fraction
    if self.problem.flow.density.type == "linear" then
        -- linear density function
        function DensityFct(w)
            return p_w + (p_s - p_w) * w
        end

        function dwDensityFct(w)
            return (p_s - p_w)
        end

        density = LuaUserFunctionNumber("DensityFct", 1);
        density:set_deriv(0, "dwDensityFct");

    elseif self.problem.flow.density.type == "exp" then
        -- exponential density function
        function DensityFct(w)
            return p_w * (p_s / p_w)^w
        end

        function dwDensityFct(w)
            return p_w * (p_s / p_w)^w * math.log(p_s / p_w)
        end

        density = LuaUserFunctionNumber("DensityFct", 1);
        density:set_deriv(0, "dwDensityFct");

    elseif self.problem.flow.density.type == "ideal" then
        -- ideal density function
        function DensityFct(w)
            return 1/(1/p_w + (1/p_s - 1/p_w) * w)
        end

        function dwDensityFct(w)
            return (-1.0 * p_w * p_s * (p_w - p_s)) / (p_s * (p_w - 1) - p_w * w)^2
        end

        density = LuaUserFunctionNumber("DensityFct", 1);
        density:set_deriv(0, "dwDensityFct");

    elseif self.problem.flow.density.type == "const" then
        density = self.problem.flow.density.min
    end

    return density
end

<<<<<<< HEAD
=======
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


>>>>>>> feature-arne
function ProblemDisc:viscosity()
    local viscosity = nil
    local mu0 = self.problem.flow.viscosity.mu0
    if self.problem.flow.viscosity.type == "const" then
        viscosity = mu0

    elseif self.problem.flow.viscosity.type == "real" then
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
