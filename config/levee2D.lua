local levee2D = 
{ 
  -- The domain specific setup
  domain = 
  {
    dim = 2,
    grid = "grids/levee2D.ugx",
    numRefs = ARGS.numRefs,
    numPreRefs = ARGS.numPreRefs,
    neededSubsets = {}
  },

  -- list of non-linear models => translated to functions
  parameter = {  
  
    { uid = "@Silt", 
      type = "vanGenuchten",
      thetaS = 0.396, thetaR = 0.131,
    --thetaS = 0.396, thetaR = 0.396,
      alpha = 0.423, n = 2.06, Ksat= 4.96e-1 }, 
    
    { uid = "@Clay",  -- modified n
      type = "vanGenuchten",
      alpha = 0.152, n = 3.06,  
      thetaS = 0.446, thetaR = 0.1, 
      Ksat= 8.2e-4 * 1e-3,},  --KSat= kappa/mu*rh0*g   <=> kappa = Ksat*mu/(rho*g)
      
    { uid = "@UserPorosity", 
      type = "const", 
      value = 1.0 }, 
      
   { uid = "@CharacteristicTime", 
      type = "const", 
      value = 1.0 }, 
      
   { uid = "@RiseTime", 
      type = "const", 
      value = 1.0 },
      
  },
  
  -- TODO: Translate table 
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
      min = 1000,       -- [ kg m^{-3} ] water density
      max = 1200        -- [ kg m^{-3} ] saltwater density
    },  
    
    viscosity = 
    { type = "real",      -- viscosity function ["const", "real"] 
      mu0 = 1e-3        -- [ kg m^{-3} ]  
    },
    air_pressure = 1013.25e2,
  },
   medium = 
   {
      {   subsets = {"CLAY"}, 
          porosity = 1.0,
          saturation = 
          { type = "vanGenuchten",
            value = "@Silt"
          },
          conductivity =
          { type  = "vanGenuchten", 
            value   = "@Silt"
          },
          diffusion   = 3.565e-6,   -- constant
          permeability  = 4.845e-13,  -- constant
      },
      {   subsets = {"SAND_LEFT","SAND_RIGHT"}, 
          porosity = 1.0,
          saturation    = 
          { type = "vanGenuchten",
            value = "@Silt"
          },
          conductivity  =
          { type      = "vanGenuchten",
            value = "@Silt"
          },
          diffusion   = 3.565e-6,
          permeability  = 4.845e-13,
      },
  },

   initial_conditions = 
   {
       { cmp = "p", value = "Levee2D_HydrostaticHead_p" },
   },

  boundary_conditions = 
  {
     {cmp = "p", type = "dirichlet", bnd = "WaterBnd", value = "Levee2D_RisingFlood_p"},
     {cmp = "p", type = "dirichlet", bnd = "ToeBnd", value = 0.0 },
     {cmp = "c", type = "dirichlet", bnd = "WaterBnd", value = 1.0},
     {cmp = "c", type = "dirichlet", bnd = "ToeBnd", value = 0.0 },
       --  {cmp = "h", type = "outflow", bnd = "AirBnd" },

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




-- Some functions (sadly global...)

Levee2D_tRise = levee2D.var.RiseTime

-- Anstieg bis auf 5.85 m
function Levee2D_Pegel(t) 
  return math.min(t/Levee2D_tRise, 1.0)*5.85 
end  

-- This is a rising flood.
function Levee2D_RisingFlood(x, z, t, si) 
  local zPegel = Levee2D_Pegel(t)
  if (z <zPegel) then return true, math.max((zPegel-z), 0.0) end
  return false, 0.0
end

function Levee2D_RisingFlood_p(x, z, t, si) 
  local t, h = Levee2D_RisingFlood(x, z, t, si)
  return t, h * 1200 * 9.81
end

-- This is an empty levee.
function Levee2D_HydrostaticHead(x, z, t, si) 
  return -z
end

function Levee2D_HydrostaticHead_p(x, z, t, si) 
  return -z * 1200 * 9.81
end


return levee2D

