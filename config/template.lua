template = 
{ 
    -- The domain specific setup
    domain = 
    {
        dim = {__type = "number", __min = 1, __max = 3}, -- problem dimension
        grid = {__type = "string"},                      -- grid name
        numRefs = {__type = "number"},                   -- refinements after distribution
        numPreRefs = {__type = "number"},                -- refinements before distribution
    },

    -- The density-driven-flow setup
    flow = 
    {
        cmp = {__type = "table", __content_type = "string"},

        gravity = {__type = "number"},    -- [ m s^{-2}�] ("standard", "no" or numeric value) 
        density =                   
        {   type = {__type = {"linear", "exp", "ideal"}},        -- density function ["linear", "exp", "ideal"]
            min = {__type = "number"},         -- [ kg m^{-3} ] water density
            max = {__type = "number"},         -- [ kg m^{-3} ] saltwater density
        },
        viscosity = 
        {   type = {__type = {"const", "real"}}, -- viscosity function ["const", "real"] 
            mu0 = {__type = "number"},           -- [ kg m^{-3} ]    
        },
    },
    medium = 
    {   __repeatable = true,
        __values = 
        {   
            {   subsets = {__type = "table", __content_type = "string"}, 
                porosity = {__type = "number"},
                saturation = 
                {   __deftypes = true,
                    __type = {"const", "exp", "vanGenuchten"},
                    -- saturation function ["const", "exp", "richards"]
                    -- if "const" define value = number
                    -- if "exp" define alpha, air_pressure and K_s
                    -- if "vanGenuchten" define a value, that fits to the uid of a parameter definition
                    __const = {
                        value = {__type = "number"},
                    },
                    __exp = {
                        value = {__type = "string"},
                    },
                    __vanGenuchten = {
                        value = {__type = "string"},
                    }
                },
                conductivity =
                {   __deftypes = true,
                    -- type         = "exp",    -- ["const", "exp", "vanGenuchten"]
                    -- if "const" define value = number
                    -- if "exp" define alpha, air_pressure, thetaS and thetaR
                    -- if "vanGenuchten" define a parameter list then value = uid
                    __type = {"const", "exp", "vanGenuchten"},
                    __const = {__type = "number"},
                    __exp = 
                    {
                        value = {__type = "string"}
                    },
                    __vanGenuchten = 
                    {
                        value = {__type = "string"},
                    },
                },
                diffusion       = {__type = "number"},
                permeability    = {__type = "number"},
            },
        },
    },
    parameter = 
    {   __repeatable = true,
        __values = 
        {
            {   uid = {__type = "string"}, 
                type = {__type = {"vanGenuchten", "exp"}},
                thetaS = {__type = "number"}, 
                thetaR = {__type = "number"},
                alpha = {__type = "number"}, 
                n = {__type = "number"}, 
                -- optional
                -- Ksat= {__type = "number"},
            },
            {   uid = {__type = "string"},
                type = {__type = {"const"}},
                value = {__type = "number"},
            },
        },
    },
    initial_conditions = 
    {   __repeatable = true,
        __values = 
        {
            {   cmp = {__type = "string"}, 
                value = {__type = {"number", "string"}} 
            },
        },
    },

    boundary_conditions = 
    {   __repeatable = true,
        __values = 
        {
            {   cmp = {__type = "string"}, 
                type = {__type = {"dirichlet"}}, 
                bnd = {__type = "string"},
                value = {__type = {"number", "string"}},
            },
        },
    },

    solverConfig = {
        solver = {
            __deftypes = true,
            __type = {"newton", "limex"},
            __newton = {
                maxSteps = {__type = "number"},
                minDef = {__type = "number"},
                reduction = {__type = "number"},
                -- optional are: 
                -- verbose = {__type = "bool"}
                -- LineSearchVerbose = {__type = "bool"}
                -- both are normaly false
            },
        },
        linSolver = {
            __deftypes = true,
            __type = {"standard", "gmg", "bicgstab"},
            __standard = {},
            __gmg = {
                maxSteps = {__type = "number"},
                minDef = {__type = "number"},
                reduction = {__type = "number"},
                },
            __bicgstab = {
                maxSteps = {__type = "number"},
                minDef = {__type = "number"},
                reduction = {__type = "number"},
                },
            }
    
  },
}