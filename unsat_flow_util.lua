-- author: Niklas Conen

-- utility functions for unsaturated density flow equations
-- problems are defined in the config file

-- namespace creation
util.unsat = util.unsat or {}


function util.unsat.conductivity(condDesc)
    local conductivity = nil
    if condDesc == "exp" then
        function ConductivityFct(w)
            return 1
        end

        function dpConductivityFct(w)
            return 1
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
        -- p = p_s + (p_w - p_s) * w
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
        -- p = p_s * (p_w / p_s)^w
        function DensityFct(w)
            return p_s * (p_w / p_s)^w
        end

        function dwDensityFct(w)
            return p_s * (p_w / p_s)^w * math.log(p_w / p_s)
        end

        density = LuaUserFunctionNumber("DensityFct", 1);
        density:set_deriv(0, "dwDensityFct");

    elseif densDesc.type == "ideal" then
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
    if satDesc.type == "const" then
        -- constant saturation
        saturation = satDesc.sat
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



-- Ansatzraum
-- Creates a approximation space for a problem
function util.unsat.CreateApproxSpace(problem, dom)
	-- the approximation space is created by using the ug4 functions
	-- documentation: line 277 in /ug4/ugcore/ugbase/lib_disc/function_spaces/approximation_space.h
	local approxSpace = ApproximationSpace(dom)
	-- adding a linear Ansatzfunktion for the variable defined in the problem file
	-- problem.flow.cmp[1] is the unknown variable
	approxSpace:add_fct("p", "Lagrange", 1)
	approxSpace:add_fct("m", "Lagrange", 1)
	approxSpace:init_levels()
	approxSpace:init_top_surface()
	approxSpace:print_statistic()
  return approxSpace
end

-- Element discretisation
function util.unsat.CreateElemDisc(subdom, dim, densDesc, poroDesc, gravDesc, satDesc, visDesc, permeability)
    -- Creates the elememt discretisation for a given medium
    -- subdom: subset with a isotropic medium
    -- dim: problem dimension
    -- density: density description
    -- porosity: medium porosity
    -- gravity: gravitational force value. force always directs downwards

	local elemDisc = {}
    -- flow equation
	elemDisc["flow"] = ConvectionDiffusion("f", subdom, "fv1")
	-- transport equation
    elemDisc["salt"] = ConvectionDiffusion("s", subdom, "fv1")

    -- constructing gravity vector
    gravity = ConstUserVector(0.0)
    gravity:set_entry(dim-1, gravDesc)

    viscosity = util.unsat.viscosity(visDesc)

    density = util.unsat.density(densDesc)

    -- hydrCond = permeability * density * gravity / viscosity

    -- hydraulic conductivity in saturated medium
    -- K_s = k * rho * g / mu
    -- ScaleAddLinkerMatrix?



    -- the permeability is the product of the hydraulic conductivity of
    -- the saturated medium K_s and the pressure dependend relative
    -- hydraulic conductivity k(p)



    local DarcyVelocity = DarcyVelocityLinker()
    DarcyVelocity:set_permeability(0)
    DarcyVelocity:set_viscosity(viscosity)
    DarcyVelocity:set_pressure_gradient(elemDisc["flow"]:gradient())
    DarcyVelocity:set_density(density)
    DarcyVelocity:set_gravity(gravity)
	print("Darcy Velocity created.")

	local storage = porosity*saturation*density

	-----------------------------------------
	-- Equation [1]
	-----------------------------------------
	-- flow equation
	-- $\partial_t (\Phi \rho_w S_w) 
	--		+ \nabla \cdot [\rho_w \vec{v}_w] = \rho_w \Gamma_w$

	local porosity = 0.5
	local saturation = 1
	local density = 1

	-- fluid storage: \Phi \rho_w S_w
	local fluidStorage = ScaleAddLinkerNumber()
	fluidStorage:add(storage)

	-- flux \rho_w \vec{v}_w
	local fluidFlux = ScaleAddLinkerVector()
	fluidFlux:add(density, DarcyVelocity)


	elemDisc["flow"]:set_mass(fluidStorage)
	elemDisc["flow"]:set_flux(-1.0*fluidFlux)
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
	diffusion:add(wDensity, -1.0)

	saturation = LuaUserFunctionNumber("exp_saturation", 1)
	saturation:set_deriv(0, "der_exp_saturation")


	saturation:set_input(elemDisc["flow"]:value())
	conductivity:set_input(elemDisc["flow"]:value())

	elemDisc["salt"]:set_mass(0.0)
	elemDisc["salt"]:set_mass_scale(fluidStorage)
	elemDisc["salt"]:set_flux(-1.0*fluidFlux)
  	elemDisc["salt"]:set_diffusion(diffusion)
end


function util.unsat.CreateDomainDisc(approxSpace, elemDisc)
	domainDisc = DomainDiscretization(approxSpace)
	domainDisc:add(elemDisc["flow"])
	domainDisc:add(elemDisc["salt"])

	--TODO: Boundaries
end

