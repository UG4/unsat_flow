-- author: Niklas Conen

-- utility functions for unsaturated density flow equations
-- problems are defined using config files

ug_load_script("ug_util.lua")

-- namespace creation
util.unsat = util.unsat or {}

local json = require("json")

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
    -- Creates the elememt discretisation for a given medium
	local elemDisc = {}
    -- flow equation
	elemDisc["flow"] = ConvectionDiffusion(self.cmp[1], subdom, "fv1")
	-- transport equation
    elemDisc["transport"] = ConvectionDiffusion(self.cmp[2], subdom, "fv1")

    local viscosity = self.mu

    local density = self.rho

    local permeability = medium.permeability

    local porosity = medium.porosity


    -- the vanGenuchten model is calculated using the Richards Plugin
    local conductivity = ProblemDisc:conductivity(medium.conductivity.value)
    local saturation = ProblemDisc:saturation(medium.saturation.value)

    -- hydraulic conductivity in saturated medium
    -- K_s = k / mu
    -- for unsaturated flow: 
    -- multiplied by relative hydraulic conductivity k(p)
    local hydrCond = ScaleAddLinkerMatrix()
    hydrCond:add(conductivity, permeability)

    -- Darcy Velocity
    -- $\vec q := -k*k(p)/mu (\grad p - \rho \vec g)$
    local DarcyVelocity = DarcyVelocityLinker()
    DarcyVelocity:set_permeability(hydrCond)
    DarcyVelocity:set_viscosity(viscosity)
    DarcyVelocity:set_pressure_gradient(self.grad_p)
    DarcyVelocity:set_density(density)
    DarcyVelocity:set_gravity(self.gravity)
	-- print("Darcy Velocity created.")

    -- molecular diffusion
    local diffusion = ScaleAddLinkerMatrix()
    diffusion:add(porosity, medium.diffusion)

	-----------------------------------------
	-- Equation [1]
	-----------------------------------------
	-- flow equation
	-- $\partial_t (\Phi \rho_w S_w) 
	--		+ \nabla \cdot [\rho_w \vec{v}_w] = \rho_w \Gamma_w$

	-- fluid storage: \Phi \rho_w S_w
    local storage = ScaleAddLinkerNumber()
	-- flux \rho_w \vec{v}_w
    local fluidFlux = ScaleAddLinkerVector()

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
        storage:add(porosity*density, saturation)
        fluidFlux:add(density, DarcyVelocity)

        elemDisc["flow"]:set_mass(storage)
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

    elemDisc["transport"]:set_mass_scale(storage)
	elemDisc["transport"]:set_velocity(fluidFlux)
	elemDisc["transport"]:set_mass(0.0)
    elemDisc["transport"]:set_diffusion(diffusion)


    -- setting inputs

    if medium.conductivity.type == "exp" then
        conductivity:set_input(0, self.capillary)

    elseif medium.conductivity.type == "vanGenuchten" then
        conductivity:set_capillary(self.capillary)
    end

    if medium.saturation.type == "exp" then
        saturation:set_input(0, self.capillary)

    elseif medium.saturation.type == "vanGenuchten" then
        saturation:set_capillary(self.capillary)
    end

    print("Created Element Discretisation for Subset ", subdom)


    -- Preparations for IO.
    local si = self.domain:subset_handler():get_subset_index(subdom)
    self.CompositeConductivity:add(si, conductivity)
    self.CompositeSaturation:add(si, saturation)
    self.CompositeDarcyVelocity:add(si, DarcyVelocity)
    self.CompositeFlux:add(si, fluidFlux)

    return elemDisc
end


function ProblemDisc:CreateDomainDisc(approxSpace)
    self.domainDisc = DomainDiscretization(approxSpace)

    local conc = GridFunctionNumberData(self.u, "c")

    -- density grid function
    self.rho = self:density()
    self.rho:set_input(0, conc)

    self.mu = self:viscosity()
    if self.problem.flow.viscosity.type ~= "const" then
        self.mu:set_input(0, conc)
    end

    self.grad_p = GridFunctionGradientData(self.u, "p")

    self.capillary = ScaleAddLinkerNumber()
    self.capillary:add(-1.0, GridFunctionNumberData(self.u, "p"))

    if self.problem.boussinesq then
        print("using boussineq approximation")
    end

    self.CompositeConductivity = CompositeUserNumber(false)
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
            self.vtk:select_element(self.rho, v)
        -- viscosity
        elseif v == "mu" and self.problem.flow.viscosity.type ~= "const" then
            self.vtk:select_element(self.mu, v)
        -- conductivity
        elseif v == "k" then
            self.vtk:select_element(self.CompositeConductivity, v)
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
end


function ProblemDisc:CreateModelMap(paramDesc)
    local modelMap = {}
    for i, medium in ipairs(paramDesc) do
        if medium.type == "vanGenuchten" then
            modelMap[medium.uid] = CreateVanGenuchtenModel(json.encode(medium))
        elseif medium.type == "const" then
            modelMap[medium.uid] = medium.value
        elseif medium.type == "exp" then
            modelMap[medium.uid] = medium
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

function ProblemDisc:conductivity(condID)
    local conductivity = nil
    local values = self.modelMap[condID]

    if type(values) == "userdata" then
        conductivity = RichardsConductivity(self.modelMap[condID])

    elseif type(values) == "table" then
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
        end
    end
    return conductivity
end


function ProblemDisc:saturation(satID)
    local saturation = nil
    local values = self.modelMap[satID]

    if type(values) == "userdata" then
        saturation = RichardsSaturation(self.modelMap[satID])
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
