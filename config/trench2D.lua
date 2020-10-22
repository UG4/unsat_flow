-- config for modelling a drainage trench with constant groundwater flow

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
    { uid = "@Silt", 
      type = "vanGenuchten",
      thetaS = 0.396, thetaR = 0.131,
      alpha = 0.423, n = 2.06, Ksat= 4.745e-6},--4.96e-1 }, 
    
    { uid = "@Clay",  -- modified n
      type = "vanGenuchten",
      alpha = 0.152, n = 3.06,  
      thetaS = 0.446, thetaR = 0.1, 
      Ksat= 8.2e-4 * 1e-3,},  --KSat= kappa/mu*rh0*g   <=> kappa = Ksat*mu/(rho*g) 
  },
  
  paramTable = { 
    ["DrainageTime"] = 1.0, 
  },

  var ={
    CharacteristicTime = 1.0,
    RiseTime = 1.0,
  },

  flow = 
  {
    type = "haline",
    cmp = {"c", "p"},

    gravity = -9.81,    -- [ m s^{-2}�] ("standard", "no" or numeric value) 
    density =           
    { type = "linear",    -- density function ["linear", "exp", "ideal"]
      min = 998.23,       -- [ kg m^{-3} ] water density
      max = 1195.0,        -- [ kg m^{-3} ] saltwater density
      w_max = 0.26,
    },  
    
    viscosity = 
    { type = "const",      -- viscosity function ["const", "real"] 
      mu0 = 1e-3        -- [ kg m^{-3} ]  
    },
    air_pressure = 1013.25e2,
  },
   medium = 
   {
      {   subsets = {"Inner"}, 
          porosity = 1.0,
          saturation = 
          { type = "vanGenuchten",
            value = "@Silt"
          },
          conductivity =
          { type  = "vanGenuchten", 
            value   = "@Silt"
          },
          diffusion   = 0.1, --3.565e-6,   -- constant
          permeability  = 4.845e-12,  -- constant
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
     {cmp = "c", type = "dirichlet", bnd = "Trench", value = 1},
     {cmp = "c", type = "dirichlet", bnd = "Aquifer", value = 0},


  },
   
  time = 
  {
    control = "limex",
    start   = 0.0,    -- [s] start time point
    stop  = 20.0,  -- [s] end time point  -- 10,000 years
    dt  = 0.1, -- [s] initial time step
    dtmin = 1e-8, -- [s] minimal time step
    dtmax = 100, -- [s] maximal time step  -- 100.0 years
    dtred = 0.5,    -- [1] reduction factor for time step
    tol   = 1e-3
  },
}


function Trench2DDrainagePressureBoundaryTime(x, y, t, tD)
  if (t <= tD) then
    return true, (2.2*t / tD - 2.0) * trench2D.flow.density.max * trench2D.flow.gravity
  else
    return true, 0.2 * trench2D.flow.density.max * trench2D.flow.gravity
  end
end

function Trench2DDrainagePressureBoundary(x, y, t)
  return Trench2DDrainagePressureBoundaryTime(x, y, t, trench2D.paramTable["DrainageTime"])  
end

function Trench2DAquiferBoundary(x, y, t)
  temp = (1.0 - y) * trench2D.flow.density.min * math.abs(trench2D.flow.gravity)
  return true, temp
end

function Trench2DPressureStart(x, y, t)
  return (1.0 - y) * trench2D.flow.density.min * math.abs(trench2D.flow.gravity)
end

return trench2D

