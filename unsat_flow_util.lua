-- author: Niklas Conen

-- utility functions for unsaturated density flow equations
-- problems are defined in the config file

-- namespace creation
util.unsat = util.unsat or {}

local json = require("json")

function util.unsat.conductivity(condDesc)
    local conductivity = nil
    local alpha = nil
    local p_a = nil
    local thetaS = nil
    local thetaR = nil

    if condDesc.type == "const" then
        conductivity = condDesc.value

    elseif condDesc.type == "exp" then
        local alpha = condDesc.alpha
        local p_a = condDesc.air_pressure
        local thetaS = condDesc.thetaS
        local thetaR = condDesc.thetaR
        function ConductivityFct(p)
            -- p is the water pressure. the pressure used by the 
            -- exponential model should be the capillary pressure (p_c).
            -- therefore the case p_a (air pressure) - p_w (water pressure) < 0
            -- needs to be handled
            p_c = p_a - p
            print(p)
            print(p_c)
            if p_c < 0 then
                return 1
            else
                return thetaR + (thetaS - thetaR) * math.exp(alpha * p)
            end
        end

        function dpConductivityFct(p)
            p_c = p_a - p
            if p_c < 0 then
                return 0
            else
                return (thetaS - thetaR) * alpha * math.exp(alpha * p)
            end
        end

        conductivity = LuaUserFunctionNumber("ConductivityFct", 1)
        conductivity:set_deriv(0, "dpConductivityFct")
    end

    return conductivity

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


function util.unsat.saturation(satDesc)
    local saturation = nil
    local alpha = nil
    local K_s = nil
    local p_a = nil

    if satDesc.type == "const" then
        -- constant saturation
        saturation = satDesc.value

    elseif satDesc.type == "exp" then
        local alpha = satDesc.alpha
        local p_a = satDesc.air_pressure

        function SaturationFct(p)
            -- p is the water pressure. the pressure used by the 
            -- exponential model should be the capillary pressure (p_c).
            -- therefore the case p_a (air pressure) - p_w (water pressure) < 0
            -- needs to be handled
            p_c = p_a - p
            if p_c < 0 then
                return 1
            else
                return K_sat * math.exp(alpha * p)
            end
        end

        function dpSaturationFct(p)
            p_c = p_a - p
            if p_c < 0 then
                return 0
            else
                return K_sat * alpha * math.exp(alpha * p)
            end
        end

        saturation = LuaUserFunctionNumber("SaturationFct", 1)
        saturation:set_deriv(0, "dpSaturationFct")
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

function util.unsat.conductivitySaturation(condDesc, satDesc, paramDesc)
    local models = {conductivity = nil, saturation = nil}
    modelMap = {}

    for i, medium in ipairs(paramDesc) do
    	if medium.type == "vanGenuchten" then
    		modelMap[medium.uid] = CreateVanGenuchtenModel(json.encode(medium))
    	end
   	end

    if condDesc.type == "vanGenuchten" then
        models["conductivity"] = RichardsConductivity(modelMap[condDesc.value])
    else
        models["conductivity"] = util.unsat.conductivity(condDesc)
    end

    if satDesc.type == "vanGenuchten" then
        models["saturation"] = RichardsSaturation(modelMap[satDesc.value])
    else
        models["saturation"] = util.unsat.saturation(condDesc)
    end

    return models
end

-- Ansatzraum
-- Creates a approximation space for a problem
function util.unsat.CreateApproxSpace(problem, dom)
	-- the approximation space is created by using the ug4 functions
	-- documentation: line 277 in /ug4/ugcore/ugbase/lib_disc/function_spaces/approximation_space.h
	local approxSpace = ApproximationSpace(dom)
	-- adding a linear Ansatzfunktion for the variable defined in the problem file
	-- problem.flow.cmp[1] is the unknown variable
	approxSpace:add_fct("p", "Lagrange", 1)
	approxSpace:add_fct("w", "Lagrange", 1)
	approxSpace:init_levels()
	approxSpace:init_top_surface()
	approxSpace:print_statistic()
  return approxSpace
end

-- Element discretisation
function util.unsat.CreateElemDisc(subdom, densDesc, porosity, gravity, satDesc, viscDesc, permDesc, condDesc, paramDesc)
    -- Creates the elememt discretisation for a given medium
	local elemDisc = {}
    -- flow equation
	elemDisc["flow"] = ConvectionDiffusion("p", subdom, "fv1")
	-- transport equation
    elemDisc["salt"] = ConvectionDiffusion("w", subdom, "fv1")

    density = util.unsat.density(densDesc)
    density:set_input(0, elemDisc["salt"]:value())

    viscosity = util.unsat.viscosity(viscDesc)
    if viscDesc.type == "real" then
        viscosity:set_input(0, elemDisc["salt"]:value())
    end

    permeability = permDesc

    -- hydraulic conductivity in saturated medium
    -- K_s = k / mu

    -- the vanGenuchten model is calculated using the Richards Plugin
    condSatModels = util.unsat.conductivitySaturation(condDesc, satDesc, paramDesc)
    conductivity = condSatModels["conductivity"]
    saturation = condSatModels["saturation"]

    if condDesc.type == "exp" then
        conductivity:set_input(0, elemDisc["flow"]:value())

    elseif condDesc.type == "vanGenuchten" then
        conductivity:set_capillary(-1.0 * elemDisc["flow"]:value())
    end

    if satDesc.type == "exp" then
        saturation:set_input(0, elemDisc["flow"]:value())

    elseif satDesc.type == "vanGenuchten" then
        saturation:set_capillary(-1.0 * elemDisc["flow"]:value())
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
    DarcyVelocity:set_gravity(gravity)
	-- print("Darcy Velocity created.")

	-----------------------------------------
	-- Equation [1]
	-----------------------------------------
	-- flow equation
	-- $\partial_t (\Phi \rho_w S_w) 
	--		+ \nabla \cdot [\rho_w \vec{v}_w] = \rho_w \Gamma_w$

	-- fluid storage: \Phi \rho_w S_w
	local storage = ScaleAddLinkerNumber()
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
	-- salt convection and diffusion => transport equation
	-- $\partial_t (\Phi \rho_w S_w \omega) 
	-- 		+ \nabla \cdot [\rho_w \omega \vec{v}_w - \rho_w D \nabla \omega] 
	--			= \rho_w \omega \Gamma_w$

	diffusion = ScaleAddLinkerMatrix()
	diffusion:add(density, 1.0)

	elemDisc["salt"]:set_mass(0.0)
	elemDisc["salt"]:set_mass_scale(storage)
	elemDisc["salt"]:set_velocity(fluidFlux)
  	elemDisc["salt"]:set_diffusion(diffusion)

    print("Created Element Discretisation for Subset ", subdom)

    return elemDisc
end


function util.unsat.CreateDomainDisc(problem, approxSpace)
	domainDisc = DomainDiscretization(approxSpace)

    local myCompositeGauge = CompositeUserNumber(false)
    local myCompositeFlux = CompositeUserVector(false)
    local myCompositeStorage = CompositeUserNumber(false)
    local myCompositeConductivity = CompositeUserNumber(false)

    -- constructing gravity vector
    local dim = problem.domain.dim
    gravity = ConstUserVector(0.0)
    gravity:set_entry(dim-1, problem.flow.gravity)

    for i,medium in ipairs(problem.medium) do -- for all media
        local mySubsetList = medium.subsets
        local elemDisc = nil

        for j, subset in ipairs(medium.subsets) do
            elemDisc = util.unsat.CreateElemDisc(   subset,
                                                    problem.flow.density,
                                                    medium.porosity,
                                                    gravity,
                                                    medium.saturation,
                                                    problem.flow.viscosity,
                                                    medium.permeability,
                                                    medium.conductivity,
                                                    problem.parameter)
            domainDisc:add(elemDisc["flow"])
            domainDisc:add(elemDisc["salt"])

            local si = domain:subset_handler():get_subset_index(mySubset)
            myCompositeGauge:add(si, elemDisc:value())
            myCompositeFlux:add(si, elemDisc:get_flux_data())
            myCompositeStorage:add(si, elemDisc:get_storage_data())
            myCompositeConductivity:add(si, elemDisc:get_conductivity_data())
        end
    end

    vtkOutput:select_element(myCompositeGauge, "gauge");
    vtkOutput:select_element(myCompositeFlux, "flux");
    vtkOutput:select_element(myCompositeStorage, "storage");
    vtkOutput:select_element(myCompositeConductivity, "conductivity");

    -- Create Boundary Conditions
    dirichletBnd = DirichletBoundary()
    for i, v in ipairs(problem.boundary_conditions) do
        if v.type == "dirichlet" then
            dirichletBnd:add(v.value, v.cmp, v.bnd)
        end
    end

    domainDisc:add(dirichletBnd)


    print("Created Domain Discretisation")

    return domainDisc

	--TODO: Boundaries
end


function util.unsat.SetInitialData(problemDesc, u)
  for index,myInitial in ipairs(problemDesc.initial_conditions) do
    Interpolate(myInitial.value, u, myInitial.cmp)
  end
end

