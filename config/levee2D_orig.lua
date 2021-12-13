-- additional constants for vanGenuchten
levee2D_rho = 1
levee2D_g = -1 -- must be negative!
rhog = 1
tstop = 1000 * 86400 -- 100 days
Levee2D_tRise = 86400

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
    { uid = "@Clay",
      type = "vanGenuchten",
      thetaS = 0.446, thetaR = 0.0,
      alpha = 0.152/rhog, n = 1.17,
      Ksat = 0.00432},

    { uid = "@Sand",
      type = "vanGenuchten",
      thetaS = 0.37, thetaR = 0.043,
      alpha = 0.087/rhog, n = 1.58,
      Ksat = 0.02592},
  },

  flow =
  {
    boussinesq = false,

    gravity = levee2D_g,      -- [ m s^{-2}], must be negative!
    density =
    { type = "const",         -- density function ["linear", "exp", "ideal"]
      min = levee2D_rho,      -- [ kg m^{-3} ] water density
      max = 1,           -- [ kg m^{-3} ] saltwater density
    },

    viscosity =
    { type = "const",          -- viscosity function ["const", "real"]
      mu0 = 1          -- [ Pa s ]
    },
    diffusion   = 18.8571e-6  -- [ m^2/s ]
  },
   medium =
   {
      {   subsets = {"CLAY"},
          porosity = "@Clay",
          saturation =
          { type = "vanGenuchten",
            value = "@Clay"
          },
          conductivity =
          { type  = "vanGenuchten",
            value   = "@Clay"
          },
      },
      {   subsets = {"SAND_LEFT","SAND_RIGHT"},
          porosity = "@Sand",
          saturation    =
          { type = "vanGenuchten",
            value = "@Sand"
          },
          conductivity  =
          { type      = "vanGenuchten",
            value = "@Sand"
          },
      },
  },

   initial =
   {
      { cmp = "c", value = 0.0 },
      { cmp = "p", value = "Levee2D_HydrostaticHead" },
   },

  boundary =
   {
      {cmp = "p", type = "dirichlet", bnd = "AirBnd", value = 0.0},
      {cmp = "p", type = "dirichlet", bnd = "WaterBnd", value = "Levee2D_HydrostaticHead"},
      --{cmp = "p", type = "dirichlet", bnd = "WaterBnd", value = "Levee2D_RisingFlood_p"},
      {cmp = "p", type = "dirichlet", bnd = "ToeBnd", value = 0.0 },
      {cmp = "c", type = "dirichlet", bnd = "ToeBnd", value = 0.0 },
      --{cmp = "c", type = "dirichlet", bnd = "WaterBnd", value = "Levee2D_RisingFlood_c"},
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
    stop	= tstop,			        -- [s]  end time point
    max_time_steps = 1000,		  -- [1]	maximum number of time steps
    dt		= 43200,		          -- [s]  initial time step
    dtmin	= ARGS.dt,	          -- [s]  minimal time step
    dtmax	= 86400,	            -- [s]  maximal time step
    dtred	= 0.5,			          -- [1]  reduction factor for time step
    tol 	= 1e-2,
  },


  -- config for vtk output
  output =
  {
    file = "./", -- must be a folder!
    data = {"c", "p", "rho", "mu", "kr", "s", "q", "ff", "tf", "af", "df", "pc", "k"},
    -- scaling factor for correct time units.
    -- 1 means all units are given in seconds
    -- if units are scaled to days, then the scaling factor should be 86400
    scale = 1
  },
}

-- rising flood
function Levee2D_RisingFlood_p(x, y, t, si)
  local pegel = math.min(t/Levee2D_tRise, 1.0)*5.85
  if (y <= pegel) then
    return true, (pegel - y)
  end
  return false, 0.0
end

function Levee2D_RisingFlood_c(x, y, t, si)
  local pegel = math.min(t/Levee2D_tRise, 1.0)*5.85
  if (y <= pegel) then
    return true, 1.0
  end
  return false, 0.0
end

-- initial pressure function
function Levee2D_HydrostaticHead(x, y, t, si)
  return (6-y)
end


return levee2D
