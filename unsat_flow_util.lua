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
    self.approxSpace:add_fct("p", "Lagrange", 1)
    self.approxSpace:add_fct("c", "Lagrange", 1)
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
    -- flow equation
	elemDisc["flow"] = ConvectionDiffusion("p", subdom, "fv1")
	-- transport equation
    elemDisc["transport"] = ConvectionDiffusion("c", subdom, "fv1")

    -- local viscosity = self.mu
    local viscosity = self:viscosity()
    if self.problem.flow.viscosity.type ~= "const" then
        viscosity:set_input(0, elemDisc["transport"]:value())
    end
    self.mu = viscosity

    local density = self:density()
    density:set_input(0, elemDisc["transport"]:value())
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
    local conductivity = ProblemDisc:conductivity(medium.conductivity.value) -- k(S)
    local saturation = ProblemDisc:saturation(medium.saturation.value) -- S

    -- volumetric fraction of the fluidphase in the porous media
    volufrac = ScaleAddLinkerNumber()
    volufrac:add(porosity, saturation)

    -- hydraulic conductivity in saturated medium is given by darcys law
    -- k_f = K*rho*g / mu
    -- for unsaturated flow:
    -- multiplied by relative hydraulic conductivity K(S)
    -- k_f = K(S) * K * rho * g / mu
    -- calculating permeability
    local permeability = nil
    local conductivityLinker = ScaleAddLinkerMatrix()
    -- if K is given by the van Genuchten modell, the permeability needs to be
    -- calculated by dividing K_sat by density times gravity and multipling with
    -- the dynamic viscosity
    if type(medium.permeability) == "string" then
        Ksat = nil
        for i, param in ipairs(self.problem.parameter) do
            if param.uid == medium.permeability then
                Ksat = param.Ksat
            end
        end

        permeability = InverseLinker()
        permeability:divide(viscosity*Ksat, (-1.0)*density*self.problem.flow.gravity)

        -- helper linker, to be able to save the permeability in CompositeUserData object
        local permcond = ScaleAddLinkerNumber()
        permcond:add(permeability, conductivity)

        conductivityLinker:add(permcond, 1)

    else
        permeability = medium.permeability
        local conductivityLinker = ScaleAddLinkerMatrix()
        conductivityLinker:add(conductivity, permeability)
    end

    -- Darcy Velocity
    -- $\vec q := -k*k(p)/mu (\grad p - \rho \vec g)$
    local DarcyVelocity = DarcyVelocityLinker()
    DarcyVelocity:set_viscosity(viscosity)
    DarcyVelocity:set_permeability(conductivityLinker)
    DarcyVelocity:set_pressure_gradient(elemDisc["flow"]:gradient())
    DarcyVelocity:set_density(density)
    DarcyVelocity:set_gravity(self.gravity)


	-----------------------------------------
	-- Equation [1]
	-----------------------------------------
	-- flow equation
	-- $\partial_t (\Phi S_w \rho_w)
	--		+ \nabla \cdot (\rho_w \vec{v}_w) = \rho_w \Gamma_w$

	-- fluid storage: \Phi S_w \rho_w
    local storage = ScaleAddLinkerNumber()
    storage:add(volufrac, density)
	-- flux of the fluid phase \rho_w \vec{v}_w
    local fluidFlux = ScaleAddLinkerVector()
    fluidFlux:add(density, DarcyVelocity)

    if self.problem.flow.boussinesq then
        -- boussinesq fluid storage: \Phi \rho_w0 S_w
        local storageB = ScaleAddLinkerNumber()
        -- boussinesq flux \rho_w0 \vec{v}_w
        local fluidFluxB = ScaleAddLinkerVector()

        storageB:add(volufrac, self.problem.flow.density.min)
        fluidFluxB:add(self.problem.flow.density.min, DarcyVelocity)

        elemDisc["flow"]:set_mass(storageB)
        elemDisc["flow"]:set_flux(fluidFluxB)

    else
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
	-- 		+ \nabla \cdot \phi S_w (\rho_w \omega q - \rho_w D \nabla \omega)
	--			= \phi S_w \rho_w \Gamma_w$

    -- diffusive flux: \phi S_w \rho_w D
    local diffusion = ScaleAddLinkerMatrix()
    diffusion:add(storage, medium.diffusion)

    local transportFlux = ScaleAddLinkerVector()
    transportFlux:add(storage, DarcyVelocity)

    elemDisc["transport"]:set_mass(0.0)
    elemDisc["transport"]:set_mass_scale(storage)
    elemDisc["transport"]:set_velocity(transportFlux)
    elemDisc["transport"]:set_diffusion(diffusion)
    elemDisc["transport"]:set_upwind(myUpwind)

    -- setting inputs
    local capillary = -1.0*elemDisc["flow"]:value()

    if medium.conductivity.type == "exp" then
        conductivity:set_input(0, capillary)

    elseif medium.conductivity.type == "vanGenuchten" then
        conductivity:set_capillary(capillary)
    end

    if medium.saturation.type == "exp" then
        saturation:set_input(0, capillary)

    elseif medium.saturation.type == "vanGenuchten" then
        saturation:set_capillary(capillary)
    end

    print("Created Element Discretisation for Subset ", subdom)


    -- Preparations for IO.
    local si = self.domain:subset_handler():get_subset_index(subdom)
    self.CompositeCapillary:add(si, capillary)
    self.CompositeConductivity:add(si, conductivity)
    self.CompositeSaturation:add(si, saturation)
    self.CompositeDarcyVelocity:add(si, DarcyVelocity)
    self.CompositeFlux:add(si, fluidFlux)
    if type(permeability) ~= "number" then
        self.CompositePermeability:add(si, permeability)
    end

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
    self.CompositeFlux = CompositeUserVector(false)
    self.CompositePermeability = CompositeUserNumber(false)

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
        if v == "p" then
           self.vtk:select_nodal(GridFunctionNumberData(self.u, v)/self.problem.output.scale^2, v)
        -- concentration
        elseif v == "c" then
            self.vtk:select_nodal(GridFunctionNumberData(self.u, v), v)
        -- density
        elseif v == "rho" then
            self.vtk:select_element(self.rho, v)
        -- viscosity
        elseif v == "mu" and self.problem.flow.viscosity.type ~= "const" then
            self.vtk:select_element(self.mu/self.problem.output.scale, v)
        -- relative conductivity
        elseif v == "kr" then
            self.vtk:select_element(self.CompositeConductivity, v)
        -- saturation
        elseif v == "s" then
            self.vtk:select_element(self.CompositeSaturation, v)
        -- darcy velocity
        elseif v == "q" then
            self.vtk:select_element(self.CompositeDarcyVelocity*self.problem.output.scale, v)
        -- transport equation flux
        elseif v == "f" then
            self.vtk:select_element(self.CompositeFlux, v)
        -- capillary pressure
        elseif v == "pc" then
            self.vtk:select(self.CompositeCapillary/self.problem.output.scale^2, v)
        -- permeability
        elseif v == "k" then
            self.vtk:select(self.CompositePermeability, v)
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
    end

    return density
end

function ProblemDisc:conductivity(condID)
    local conductivity = nil
    local values = self.modelMap[condID]

    if type(values) == "userdata" then
        conductivity = RichardsConductivity(self.modelMap[condID])
    end
    return conductivity
end


function ProblemDisc:saturation(satID)
    local saturation = nil
    local values = self.modelMap[satID]

    if type(values) == "userdata" then
        saturation = RichardsSaturation(self.modelMap[satID])
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
