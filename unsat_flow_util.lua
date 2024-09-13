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
    self.gravity:set_entry(problemDesc.domain.dim - 1, problemDesc.flow.gravity)

    self.modelMap = nil
    if problemDesc.parameter ~= nil then
        print(problemDesc.parameter)
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

    elemDisc["flow"] = ConvectionDiffusion("p", subdom, "fv1")      -- flow equation
    elemDisc["transport"] = ConvectionDiffusion("c", subdom, "fv1") -- transport equation

    -- ##### Density #####
    local density = self:density()
    if self.problem.flow.density.type ~= "const" then
        density:set_input(0, elemDisc["transport"]:value())
    end
    self.rho = density

    -- ##### Porosity #####
    local porosity = medium.porosity -- phi
    -- if the porosity is given as a parameter it must be equal to the saturated water content
    if type(medium.porosity) == "string" then
        for i, param in ipairs(self.problem.parameter) do
            if param.uid == medium.porosity then
                porosity = param.thetaS -- the porosity is equal to the saturated water content
            end
        end
    end

    -- the vanGenuchten model is calculated using the Richards Plugin
    -- Capillary pressure (auxiliary).
    -- air pressure is set to 0 => p_c = - p_w
    local capillary = ScaleAddLinkerNumber()
    capillary:add(-1.0, elemDisc["flow"]:value())

    self:assert_richards_parameters(medium)
    local RichardsUserData = RichardsUserDataFactory(capillary)

    local conductivity = ConstUserNumber(1.0)
    if (type(medium.conductivity.value) == "string") then 
        conductivity = RichardsUserData:create_conductivity(self.modelMap[medium.conductivity.value])
    end

    local saturation = ConstUserNumber(1.0)
    print(type(medium.saturation.value))
    if (type(medium.saturation.value) == "number") then
        saturation = ConstUserNumber(medium.saturation.value)
    elseif (type(medium.saturation.value) == "string") then
        saturation = RichardsUserData:create_saturation(self.modelMap[medium.saturation.value])
    end

    local DarcyVelocity = DarcyVelocityLinker()
    -- to calculate the relative hydraulic conductivity one has to specify either the mediums
    -- permeability K and viscosity mu or the hydraulic conductivity in a saturated medium Ksat
    if medium.permeability and self.problem.flow.viscosity then 
        -- case 1: permeability and viscosity are given
        -- Darcy Velocity will be: $\vec q := -k_r*K/mu (\grad p - \rho \vec g)$
        local permeability = medium.permeability
        local khyd = ScaleAddLinkerMatrix()
        khyd:add(conductivity, permeability) -- Here Ksat must be 1.0!

        -- ##### Viscosity #####
        local viscosity = self:viscosity()
        if type(viscosity) ~= "number" then
            if self.problem.flow.viscosity.type ~= "const" then
                viscosity:set_input(0, elemDisc["transport"]:value())
            end
        end
        
        self.mu = viscosity

        DarcyVelocity:set_permeability(khyd)
        DarcyVelocity:set_viscosity(viscosity)
    else -- case 2: Ksat is given
        -- $\vec q := -k_r*Ksat/rho*g (\grad p - \rho \vec g)$
        local permeability = ScaleAddLinkerMatrix()
        permeability:add(conductivity, 1.0)
        DarcyVelocity:set_permeability(permeability)
        DarcyVelocity:set_viscosity(density * math.abs(self.problem.flow.gravity))
    end

    DarcyVelocity:set_pressure_gradient(elemDisc["flow"]:gradient())
    DarcyVelocity:set_density(density)
    DarcyVelocity:set_gravity(self.gravity)

    local volufrac = ScaleAddLinkerNumber() -- theta
    volufrac:add(saturation, 1.0)

    -----------------------------------------
    -- Equation [1]
    -----------------------------------------
    -- flow equation
    -- $\partial_t (\Phi S_w \rho_w)
    --		+ \nabla \cdot (\rho_w \vec{v}_w) = \rho_w \Gamma_w$

    -- fluid storage: \Phi \rho_w S_w
    local fluidStorage = ScaleAddLinkerNumber()
    -- flux \rho_w \vec{v}_w
    local fluidFlux = ScaleAddLinkerVector()

    -- fluidStorage:add(porosity*density, saturation) -- INSERT: saturation  or 1.0
    fluidStorage:add(volufrac, density)
    fluidFlux:add(density, DarcyVelocity)

    -- oberbeck-boussinesq approximation
    if self.problem.flow.boussinesq then
        -- boussinesq fluid storage: \Phi \rho_w0 S_w
        local fluidStorageB = ScaleAddLinkerNumber()
        -- boussinesq flux \rho_w0 \vec{v}_w
        local fluidFluxB = ScaleAddLinkerVector()

        fluidStorageB:add(self.problem.flow.density.min, volufrac)
        fluidFluxB:add(self.problem.flow.density.min, DarcyVelocity)

        elemDisc["flow"]:set_mass(fluidStorageB)
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
    -- 		+ \nabla \cdot (\rho_w \omega q - \rho_w,min D \nabla \omega)
    --			= \phi S_w \rho_w \Gamma_w$

    -- where  D = (\Phi S_w) Dmol + Dmech(q)

    -- molecular diffusion
    local transportDiffTensor = ScaleAddLinkerMatrix()
    transportDiffTensor:add(volufrac * density, self.problem.flow.diffusion)

    -- Scheidegger dispersion.
    if (medium.alphaT) then
        local dispersion = BearScheidegger(DarcyVelocity, medium.alphaT, medium.alphaL)
        transportDiffTensor:add(density, dispersion)
    end

    elemDisc["transport"]:set_mass(0.0)
    elemDisc["transport"]:set_mass_scale(fluidStorage) -- * w
    elemDisc["transport"]:set_velocity(fluidFlux)
    elemDisc["transport"]:set_diffusion(transportDiffTensor)
    
    if self.problem.flow.upwind == "full" then
        elemDisc["transport"]:set_upwind(FullUpwind())
    elseif self.problem.flow.upwind == "partial" then
        elemDisc["transport"]:set_upwind(PartialUpwind())
    end

    print("Created Element Discretisation for Subset ", subdom)

    -- Preparations for IO.
    local advFlux = ScaleAddLinkerVector()
    c_value = GridFunctionNumberData(self.u, "c")
    advFlux:add(c_value, fluidFlux)

    local gradC = GridFunctionGradientData(self.u, "c")
    local difFlux = ScaleAddLinkerVector()
    difFlux:add(volufrac, self.problem.flow.density.min * self.problem.flow.diffusion * gradC)

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
    self.CompositeSaltMass:add(si, volufrac * density * elemDisc["transport"]:value())

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

    for i, medium in ipairs(self.problem.medium) do -- for all media
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
        for i, mySource in ipairs(self.problem.sources) do -- for all sources
            -- Read Vector.
            local mycoord = Vec()

            mycoord:set_coord(0, mySource.coord[1])
            mycoord:set_coord(1, mySource.coord[2])
            if (dim == 3) then mycoord:set_coord(2, mySource.coord[3]) end
            -- Add primary.
            local elemDisc = DiracSourceDisc(mySource.cmp, mySource.subset)
            elemDisc:add_source(mySource.strength, mycoord)
            self.domainDisc:add(elemDisc)

            print("Added point source for " .. mySource.cmp .. " at ")
            print(mySource.coord)
            print("with strength " .. mySource.strength)

            -- Add source for substances
            for j, mySub in ipairs(mySource.substances) do -- for all components
                local elemDiscSub = DiracSourceDisc(mySub.cmp, mySource.subset)
                elemDiscSub:add_transport_sink(mySource.strength)
                self.domainDisc:add(elemDiscSub)
                print("Added point source for " .. mySub.cmp .. " at ")
                print(mySource.coord)
                print("with strength " .. mySource.strength)
            end
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

        if v.type == "flux" or v.type =="neumann" then
            -- Neumann-type
            neumannBnd[v.cmp] = neumannBnd[v.cmp] or NeumannBoundary(v.cmp, "fv1")
            neumannBnd[v.cmp]:add(v.value, v.bnd, v.inner)
            print("Added Neumann Boundary with value " .. v.value .. " for " .. v.cmp .. " on subset " .. v.bnd)
        end
    end

    if (dirichletBnd) then self.domainDisc:add(dirichletBnd) end
    if (neumannBnd["p"]) then self.domainDisc:add(neumannBnd["p"]) end
    if (neumannBnd["c"]) then self.domainDisc:add(neumannBnd["c"]) end


    local constraints = OneSideP1Constraints()
    self.domainDisc:add(constraints)

    print("Created Domain Discretisation")
    return self.domainDisc
end

function ProblemDisc:assert_richards_parameters(medium)
    -- if the porosity is given as a parameter it must be equal to the saturated water content
    local porosity = medium.porosity

    local cuid = medium.conductivity.value
    local suid = medium.saturation.value

    assert(cuid == suid, "conductivity and saturation must be defined for the same medium")

    local Ksat = self:lookup(cuid, "Ksat")
    
    local case1 = (Ksat == nil) and (medium.permeability ~= nil and self.problem.flow.viscosity ~= nil)
    local case2 = (Ksat ~= nil) and (medium.permeability == nil and self.problem.flow.viscosity == nil)
    print("Permeability and Viscosity are given, but not hydraulic conductivity: " .. tostring(case1))
    print("Only hydraulic conductivity is given: " .. tostring(case2))
    assert(case1 or case2, "either Ksat or permeability and viscosity must be defined")
end

function ProblemDisc:lookup(uid, param)
    -- looks up values by key
    local value = nil
    for i, p in ipairs(self.problem.parameter) do
        if p.uid == uid then
            value = p[param]
        end
    end
    return value
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
        elseif (v == "mu" or v == "all") and self.problem.flow.viscosity then
            if type(self.problem.flow.viscosity) ~= "number" and self.problem.flow.viscosity.type ~= "const" then
                self.vtk:select_element(self.mu, v)
            end
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

-- Create a map for "string" -> model
function ProblemDisc:CreateModelMap(paramDesc)
    local modelMap = {}
    local mfactory = RichardsModelFactory()

    for i, medium in ipairs(paramDesc) do
        if medium.type == "vanGenuchten" then
            modelMap[medium.uid] = mfactory:create_van_genuchten(json.encode(medium))
            print(medium.uid .. ":v->" .. modelMap[medium.uid]:config_string())
        elseif medium.type == "exp" then
            modelMap[medium.uid] = mfactory:create_exponential(json.encode(medium))
            print(medium.uid .. ":e->" .. modelMap[medium.uid]:config_string())
        elseif medium.type == "const" then
            modelMap[medium.uid] = medium.value
            print(medium.uid .. ":c->" .. modelMap[medium.uid])
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
            return p_w * (p_s / p_w) ^ w
        end

        function dwDensityFct(w)
            return p_w * (p_s / p_w) ^ w * math.log(p_s / p_w)
        end

        density = LuaUserFunctionNumber("DensityFct", 1);
        density:set_deriv(0, "dwDensityFct");
    elseif self.problem.flow.density.type == "ideal" then
        -- ideal density function
        function DensityFct(w)
            return 1 / (1 / p_w + (1 / p_s - 1 / p_w) * w)
        end

        function dwDensityFct(w)
            return (-1.0 * p_w * p_s * (p_w - p_s)) / (p_s * (p_w - 1) - p_w * w) ^ 2
        end

        density = LuaUserFunctionNumber("DensityFct", 1);
        density:set_deriv(0, "dwDensityFct");
    elseif self.problem.flow.density.type == "const" then
        density = self.problem.flow.density.min
    end

    return density
end

function ProblemDisc:viscosity()
    local viscosity = nil
    if type(self.problem.flow.viscosity) == "number" then
        return self.problem.flow.viscosity
    end

    local mu0 = self.problem.flow.viscosity.mu0
    if self.problem.flow.viscosity.type == "const" then
        viscosity = mu0
    elseif self.problem.flow.viscosity.type == "real" then
        function ViscosityFct(w)
            return mu0 * (1 + 1.85 * w - 4.1 * w ^ 2 + 44.5 * w ^ 3)
        end

        function dwViscosityFct(w)
            return mu0 * (1.85 - 8.2 * w + 133.5 * w ^ 2)
        end

        viscosity = LuaUserFunctionNumber("ViscosityFct", 1);
        viscosity:set_deriv(0, "dwViscosityFct");
    end

    return viscosity
end
