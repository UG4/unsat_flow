-- additional constants for vanGenuchten
levee2D_rho = 998.23
levee2D_g = -7.323e10 -- must be negative!
rhog = (-1.0)*levee2D_rho*levee2D_g
tstop = 100 -- 100 days
Levee2D_tRise = 1

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
      thetaS = 0.250, thetaR = 0.153,
      alpha = 0.79/rhog, n = 10.4,
      Ksat = 1.08},

    { uid = "@TouchetSiltLoam",
      type = "vanGenuchten",
      thetaS = 0.469, thetaR = 0.190,
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

  flow =
  {
    boussinesq = false,

    gravity = levee2D_g,  -- [ m day^{-2}], must be negative!
    density =
    { type = "ideal",     -- density function ["linear", "exp", "ideal"]
      min = levee2D_rho,  -- [ kg m^{-3} ] water density
      max = 1025.0,       -- [ kg m^{-3} ] saltwater density
    },

    viscosity =
    { type = "real",      -- viscosity function ["const", "real"]
      mu0 = 1.16e-8       -- [ Pa d ]
    },
    diffusion   = 1.63    -- [ m^2/day ]
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
    {cmp = "p", type = "dirichlet", bnd = "WaterBnd", value = "Levee2D_RisingFlood_p"},
    {cmp = "p", type = "dirichlet", bnd = "ToeBnd", value = 0.0 },
    {cmp = "c", type = "dirichlet", bnd = "ToeBnd", value = 0.0 },
    {cmp = "c", type = "dirichlet", bnd = "WaterBnd", value = "Levee2D_RisingFlood_c"},
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
      dt		= 0.5,		            -- [s]  initial time step
      dtmin	= ARGS.dt,	          -- [s]  minimal time step
      dtmax	= 1.0,	              -- [s]  maximal time step
      dtred	= 0.5,			          -- [1]  reduction factor for time step
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
    scale = 86400
  },
}

-- rising flood
function Levee2D_RisingFlood_p(x, y, t, si)
  local pegel = math.min(t/Levee2D_tRise, 1.0)*5.85
  if (y <= pegel) then
    return true, (pegel - y) * rhog
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
  return (6-y)*rhog
end


return levee2D