-- config for modelling a drainage trench with constant groundwater flow

trench2Drho = 1.0
trench2Dg = -9.81
trench2Drhog = -1.0 * trench2Drho * trench2Dg

local trench2D = 
{ 
  -- The domain specific setup
  domain = 
  {
    dim = 2,
    grid = "grids/trench2D.ugx",
    numRefs = ARGS.numRefs,
    numPreRefs = ARGS.numPreRefs,
  },

  -- list of non-linear models => translated to functions
  parameter = {  
    { uid = "@SiltLoamP", 
      type = "vanGenuchten",
      thetaS = 0.396/0.396, thetaR = 0.131/0.396,
      alpha = 0.423/(trench2Drhog), n = 2.06 },
  },
  
  paramTable = { 
    ["DrainageTime"] = 0.2, 
  },

  flow = 
  {
    type = "haline",
    cmp = {"p", "c"},

    gravity = trench2Dg,    -- [ m s^{-2}�] ("standard", "no" or numeric value) 
    density =           
    { type = "linear",    -- density function ["linear", "exp", "ideal"]
      min = 1.0,       -- [ kg m^{-3} ] water density
      max = 1.2,        -- [ kg m^{-3} ] saltwater density
      w_max = 0.26,
    },  
    
    viscosity = 
    { type = "const",      -- viscosity function ["const", "real"] 
      mu0 = 1e-3        -- [ kg m^{-3} ]  
    },
    air_pressure = 1,
  },
   medium = 
   {
      {   subsets = {"Inner"}, 
          porosity = 1.0,
          saturation = 
          { type = "vanGenuchten",
            value = "@SiltLoamP"
          },
          conductivity =
          { type  = "vanGenuchten", 
            value   = "@SiltLoamP"
          },
          diffusion  = 0.1,--3.565e-6 * 0.396,   -- constant
          permeability  = (0.0496 * 1e-3) / (-9.81*trench2Drho),  -- constant
      },
  },

   initial_conditions = 
   {
       { cmp = "p", value = "Trench2DPressureStart" },
       { cmp = "c", value = 0}
   },

  boundary_conditions = 
  {
     {cmp = "p", type = "dirichlet", bnd = "Trench", value = "Trench2DDrainagePressureBoundary"},
     {cmp = "p", type = "dirichlet", bnd = "Aquifer", value = "Trench2DAquiferBoundary" },
     {cmp = "c", type = "dirichlet", bnd = "Trench", value = 1.0},
     {cmp = "c", type = "dirichlet", bnd = "Aquifer", value = 0},


  },

  solverConfig = {
    solver = {
      type = "newton",
      maxSteps = 10,
      minDef = 5e-8,
      reduction = 1e-10,
      verbose = false,
      LineSearchVerbose = false,
    },
    linSolver = {
      type = "standard",
      maxSteps = 60,
      minDef = 1e-8,
      reduction = 1e-8,
    },
  },

   
  time = 
  {
    control = "limex",
    start   = 0.0,    -- [s] start time point
    stop  = 10.0,  -- [s] end time point  -- 10,000 years
    dt  = 0.01, -- [s] initial time step
    dtmin = 1e-8, -- [s] minimal time step
    dtmax = 100, -- [s] maximal time step  -- 100.0 years
    dtred = 0.5,    -- [1] reduction factor for time step
    tol   = 1e-3
  },
}


function Trench2DDrainagePressureBoundaryTime(x, y, t, tD)
  if (t <= tD) then
    return true, (2.2*t / tD - 2.0)
  else
    return true, 0.2
  end
end

function Trench2DDrainagePressureBoundary(x, y, t)
  b, p = Trench2DDrainagePressureBoundaryTime(x, y, t, trench2D.paramTable["DrainageTime"])  
  return b, p * trench2Drhog
end

function Trench2DAquiferBoundary(x, y, t)
  return true, (1.0 - y) * trench2Drhog
end

function Trench2DPressureStart(x, y, t)
  return (1.0 - y) * trench2Drhog
end

return trench2D

