-- additional constants for vanGenuchten
levee2D_rho = 998.23
levee2D_g = -9.81 -- must be negative!
rhog = (-1.0)*levee2D_rho*levee2D_g

local levee2D =
{
  -- The domain specific setup
  domain =
  {
    dim = 2,
    grid = "grids/levee2D.ugx",
    numRefs = ARGS.numRefs,
    numPreRefs = ARGS.numPreRefs,
  },

  -- medium parameters for vanGenuchten Model
  parameter = {
    { uid = "@Sandstone",
      type = "vanGenuchten",
      thetaS = 0.153, thetaR = 0.250,
      alpha = 0.79/rhog, n = 10.4,
      Ksat = 1.08},

    { uid = "@TouchetSiltLoam",
      type = "vanGenuchten",
      thetaS = 0.190, thetaR = 0.469,
      alpha = 0.50/rhog, n = 7.09,
      Ksat = 3.03},

    { uid = "@SiltLoam",
      type = "vanGenuchten",
      thetaS = 0.396, thetaR = 0.131,
      alpha = 0.423/rhog, n = 2.06,
      Ksat = 0.0496},

    { uid = "@Clay",
      type = "vanGenuchten",
      thetaS = 0.446, thetaR = 0.0,
      alpha = 0.152/rhog, n = 1.17,
      Ksat = 8.2e-4}
    },

  },

  var ={
    CharacteristicTime = 1.0,
    RiseTime = 1.0,
  },

  flow =
  {
    boussinesq = false,

    gravity = levee2D_g,      -- [ m s^{-2}], must be negative!
    density =
    { type = "ideal",         -- density function ["linear", "exp", "ideal"]
      min = levee2D_rho,      -- [ kg m^{-3} ] water density
      max = 1025.0,           -- [ kg m^{-3} ] saltwater density
    },

    viscosity =
    { type = "real",          -- viscosity function ["const", "real"]
      mu0 = 1.002e-3          -- [ kg m^{-3} ]
    },
    diffusion   = 18.8571e-6  -- constant
  },
   medium =
   {
      {   subsets = {"CLAY"},
          porosity = "@SiltLoam",
          saturation =
          { type = "vanGenuchten",
            value = "@SiltLoam"
          },
          conductivity =
          { type  = "vanGenuchten",
            value   = "@SiltLoam"
          },
          permeability  = "@SiltLoam",    -- uid of a medium defined under parameter or number
      },
      {   subsets = {"SAND_LEFT","SAND_RIGHT"},
          porosity = "@SiltLoam",
          saturation    =
          { type = "vanGenuchten",
            value = "@SiltLoam"
          },
          conductivity  =
          { type      = "vanGenuchten",
            value = "@SiltLoam"
          },
          permeability  = "@SiltLoam",    -- uid of a medium defined under parameter or number
      },
  },

   initial =
   {
      { cmp = "c", value = 0.0 },
      { cmp = "p", value = "Levee2D_HydrostaticHead_p" },
   },

  boundary =
  {
     {cmp = "p", type = "dirichlet", bnd = "RightBnd", value = 0},
     {cmp = "p", type = "dirichlet", bnd = "WaterBnd", value = "Levee2D_RisingFlood_p"},
     {cmp = "p", type = "dirichlet", bnd = "ToeBnd", value = 0.0 },
     {cmp = "c", type = "dirichlet", bnd = "WaterBnd", value = 1.0},
     {cmp = "c", type = "dirichlet", bnd = "ToeBnd", value = 0.0 },
  },

  linSolver =
  { type = "bicgstab",			                      -- linear solver type ["bicgstab", "cg", "linear"]
    precond =
    { type 		= "gmg",	                          -- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
      smoother 	= {type = "ilu", overlap = true},	-- gmg-smoother ["ilu", "ilut", "jac", "gs", "sgs"]
      cycle		= "V",		                          -- gmg-cycle ["V", "F", "W"]
      preSmooth	= 3,		                          -- number presmoothing steps
      postSmooth 	= 3,		                        -- number postsmoothing steps
      rap			= true,		                          -- comutes RAP-product instead of assembling if true
      baseLevel	= ARGS.numPreRefs,                -- gmg - baselevel
    },
    convCheck =
      { type		= "standard",
        iterations	= 30,		-- number of iterations
        absolute	= 0.5e-8,	-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
        reduction	= 1e-7,		-- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
        verbose		= true		-- print convergence rates if true
      }
  },

  time =
  {
      control	= "limex",
      start 	= 0.0,				      -- [s]  start time point
      stop	= 1000.0,			        -- [s]  end time point
      max_time_steps = 10000,		  -- [1]	maximum number of time steps
      dt		= tstop/50,		        -- [s]  initial time step
      dtmin	= ARGS.dt,	          -- [s]  minimal time step
      dtmax	= 10.0,	              -- [s]  maximal time step
      dtred	= 0.1,			          -- [1]  reduction factor for time step
      tol 	= 1e-2,
  },


  -- config for vtk output
  output =
  {
    file = "simulations/levee2D", -- must be a folder!
    data = {"c", "p", "rho", "mu", "kr", "s", "q", "ff", "tf", "af", "df", "pc", "k"},
    -- scaling factor for correct time units.
    -- 1 means all units are given in seconds
    -- if units are scaled to days, then the scaling factor should be 86400
    scale = 1
  },
}

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
