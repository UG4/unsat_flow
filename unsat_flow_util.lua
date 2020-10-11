-- author: Niklas Conen

-- utility functions for unsaturated density flow equations
-- problems are defined in the config file

-- namespace creation
util.unsat = util.unsat or {}

local json = require("json")

-- object oriented clustering of the Approximation Space, Element Discretisation and
-- Domain Discretisation
ProblemDisc = {}

function ProblemDisc:new(problemDesc, dom, vtk)
    assert(problemDesc, "No Input defined")
    problemObj = {}
    setmetatable(problemObj, self)
    self.__index = self

    self.problem = problemDesc
    self.domain = dom
    self.vtk = vtk
    self.cmp = problemDesc.flow.cmp
    self.approxSpace = nil
    self.ElemDisc = nil
    self.domainDisc = nil

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
    approxSpace = ApproximationSpace(self.domain)
    approxSpace:add_fct(self.cmp[1], "Lagrange", 1)
    approxSpace:add_fct(self.cmp[2], "Lagrange", 1)
    approxSpace:init_levels()
    approxSpace:init_top_surface()
    approxSpace:print_statistic()
    self.approxSpace = approxSpace
  return approxSpace
end

-- Element discretisation
function ProblemDisc:CreateElemDisc(subdom, medium)
    -- Creates the elememt discretisation for a given medium
	elemDisc = {}
    -- flow equation
	elemDisc["flow"] = ConvectionDiffusion(self.cmp[1], subdom, "fv1")
	-- transport equation
    elemDisc["transport"] = ConvectionDiffusion(self.cmp[2], subdom, "fv1")

    density = util.unsat.density(self.problem.flow.density)
    density:set_input(0, elemDisc["transport"]:value())

    viscosity = util.unsat.viscosity(self.problem.flow.viscosity)
    if self.problem.flow.viscosity.type == "real" then
        viscosity:set_input(0, elemDisc["transport"]:value())
    end

    permeability = medium.permeability

    porosity = medium.porosity
    -- hydraulic conductivity in saturated medium
    -- K_s = k / mu

    -- the vanGenuchten model is calculated using the Richards Plugin
    local conductivity = ProblemDisc:conductivity(medium.conductivity.value)
    local saturation = ProblemDisc:saturation(medium.saturation.value)

    capillary = self.problem.flow.air_pressure - elemDisc["flow"]:value()

    if medium.conductivity.type == "exp" then
        conductivity:set_input(0, capillary)

    elseif medium.conductivity.type == "vanGenuchten" then
        conductivity:set_capillary(-1.0 * capillary)
    end

    if medium.saturation.type == "exp" then
        saturation:set_input(0, capillary)

    elseif medium.saturation.type == "vanGenuchten" then
        saturation:set_capillary(-1.0 * capillary)
    end

    -- the permeability is multiplied by relative
    -- hydraulic conductivity k(p)
    hydrCond = ScaleAddLinkerMatrix()
    hydrCond:add(conductivity, permeability)

    -- Darcy Velocity
    -- $\vec q := -k*k(p)/mu (\grad p - \rho \vec g)$
    local DarcyVelocity = DarcyVelocityLinker()
    DarcyVelocity:set_permeability(hydrCond)
    DarcyVelocity:set_viscosity(viscosity)
    DarcyVelocity:set_pressure_gradient(elemDisc["flow"]:gradient())
    DarcyVelocity:set_density(density)
    DarcyVelocity:set_gravity(self.gravity)
	-- print("Darcy Velocity created.")

	-----------------------------------------
	-- Equation [1]
	-----------------------------------------
	-- flow equation
	-- $\partial_t (\Phi \rho_w S_w) 
	--		+ \nabla \cdot [\rho_w \vec{v}_w] = \rho_w \Gamma_w$

	-- fluid storage: \Phi \rho_w S_w
	local storage = ScaleAddLinkerNumber()
    print(porosity, density, saturation)
    storage:add(porosity*density, saturation)

	-- flux \rho_w \vec{v}_w
	local fluidFlux = ScaleAddLinkerVector()
	fluidFlux:add(density, DarcyVelocity)

	elemDisc["flow"]:set_mass(storage)
	elemDisc["flow"]:set_flux(fluidFlux)
	elemDisc["flow"]:set_mass_scale(0.0)
  	elemDisc["flow"]:set_diffusion(0.0)

	-----------------------------------------
	-- Equation [2]
	-----------------------------------------
	-- transport convection and diffusion => transport equation
	-- $\partial_t (\Phi \rho_w S_w \omega) 
	-- 		+ \nabla \cdot [\rho_w \omega \vec{v}_w - \rho_w D \nabla \omega] 
	--			= \rho_w \omega \Gamma_w$

	diffusion = ScaleAddLinkerMatrix()
	diffusion:add(density, 1.0)

	elemDisc["transport"]:set_mass(0.0)
	elemDisc["transport"]:set_mass_scale(storage)
	elemDisc["transport"]:set_velocity(fluidFlux)
  	elemDisc["transport"]:set_diffusion(diffusion)

    print("Created Element Discretisation for Subset ", subdom)

    return elemDisc
end


function ProblemDisc:CreateDomainDisc(approxSpace)
	domainDisc = DomainDiscretization(approxSpace)

    local myCompositeGauge = CompositeUserNumber(false)
    local myCompositeFlux = CompositeUserVector(false)
    local myCompositeStorage = CompositeUserNumber(false)
    local myCompositeConductivity = CompositeUserNumber(false)

    for i,medium in ipairs(self.problem.medium) do -- for all media
        local elemDisc = nil

        for j, subset in ipairs(medium.subsets) do
            elemDisc = ProblemDisc:CreateElemDisc(subset, medium)

            domainDisc:add(elemDisc["flow"])
            domainDisc:add(elemDisc["transport"])

            --local si = domain:subset_handler():get_subset_index(mySubset)
            --myCompositeGauge:add(si, elemDisc:value())
            --myCompositeFlux:add(si, elemDisc:get_flux_data())
            --myCompositeStorage:add(si, elemDisc:get_storage_data())
            --myCompositeConductivity:add(si, elemDisc:get_conductivity_data())
        end
    end

    --vtkOutput:select_element(myCompositeGauge, "gauge");
    --vtkOutput:select_element(myCompositeFlux, "flux");
    --vtkOutput:select_element(myCompositeStorage, "storage");
    --vtkOutput:select_element(myCompositeConductivity, "conductivity");

    -- Create Boundary Conditions
    dirichletBnd = DirichletBoundary()
    for i, v in ipairs(self.problem.boundary_conditions) do
        if v.type == "dirichlet" then
            dirichletBnd:add(v.value, v.cmp, v.bnd)
        end
    end

    domainDisc:add(dirichletBnd)


    print("Created Domain Discretisation")
    self.domainDisc = domainDisc
    return domainDisc
end

function ProblemDisc:CreateModelMap(paramDesc)
    modelMap = {}
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

function util.unsat.SetInitialData(problemDesc, u)
  for index,myInitial in ipairs(problemDesc.initial_conditions) do
    Interpolate(myInitial.value, u, myInitial.cmp)
  end
end

function util.unsat.density(densDesc)   
    local p_w = densDesc.min
    local p_s = densDesc.max
    local density = nil
    
    -- p_s: max density -> saline
    -- p_w: min density -> water
    -- w: mass fraction
    if densDesc.type == "linear" then
        -- linear density function
        function DensityFct(w)
            return p_s + (p_w - p_s) * w
        end

        function dwDensityFct(w)
            return p_w - p_s
        end

        density = LuaUserFunctionNumber("DensityFct", 1);
        density:set_deriv(0, "dwDensityFct");

    elseif densDesc.type == "exp" then
        -- exponential density function
        function DensityFct(w)
            return p_s * (p_w / p_s)^w
        end

        function dwDensityFct(w)
            return p_s * (p_w / p_s)^w * math.log(p_w / p_s)
        end

        density = LuaUserFunctionNumber("DensityFct", 1);
        density:set_deriv(0, "dwDensityFct");

    elseif densDesc.type == "ideal" then
        -- ideal density function
        function DensityFct(w)
            return 1/(1/p_s + (1/(p_w - p_s))^w)
        end

        function dwDensityFct(w)
            return -1 * ((1/p_w - 1/p_s)^w * math.log(1/p_w - 1/p_s))/((1/p_w - 1/p_s)^w + 1/p_s)^2
        end

        density = LuaUserFunctionNumber("DensityFct", 1);
        density:set_deriv(0, "dwDensityFct");
    end

    return density
end

function ProblemDisc:conductivity(condID)
    conductivity = nil
    values = self.modelMap[condID]

    if type(values) == "userdata" then -- not the prettiest solution. this should be reworked
        conductivity = RichardsConductivity(modelMap[condID])
    elseif type(values) == "table" then
        if values.type == "exp" then
            local alpha = values.alpha
            local thetaS = values.thetaS
            local thetaR = values.thetaR
            function ConductivityFct(p)
                if p < 0 then
                    return 1
                else
                    return thetaR + (thetaS - thetaR) * math.exp(alpha * p)
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
    saturation = nil
    values = self.modelMap[satID]

    if type(values) == "userdata" then
        saturation = RichardsSaturation(modelMap[satID])
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

function util.unsat.viscosity(visDesc)
    local viscosity = nil
    local mu0 = nil
    if visDesc.type == "const" then
        viscosity = visDesc.mu0

    elseif visDesc.type == "real" then
        mu0 = visDesc.mu0

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
