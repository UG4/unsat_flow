-- additional constants for vanGenuchten
henry2D_rho = 1000
henry2D_g = -9.8 -- must be negative!
rhog = (-1.0)*henry2D_rho*henry2D_g
tstop = 50 * 86400 -- 100 days


local henry =
{
  -- The domain specific setup
  domain =
  {
    dim = 2,
    grid = "grids/henry.ugx",
    numRefs = ARGS.numRefs,
    numPreRefs = ARGS.numPreRefs,
  },

  -- medium parameters for vanGenuchten Model
  parameter = {
    { uid = "@SiltLoam",
      type = "vanGenuchten",
      thetaS = 0.35, thetaR = 0.131,
      alpha = 0.423/rhog, n = 2.06,
      Ksat = 1e-2},
  },

  flow =
  {
    boussinesq = false,

    gravity = henry2D_g,      -- [ m s^{-2}], must be negative!
    density =
    { type = "linear",        -- density function ["linear", "exp", "ideal"]
      min = henry2D_rho,      -- [ kg m^{-3} ] water density
      max = 1025,          -- [ kg m^{-3} ] saltwater density
    },
    diffusion   = 3.77143e-6,   -- [ m^2/s ]
    upwind  = "full"
  },

  medium =
   {
      {   subsets = {"Medium"},
          porosity = "@SiltLoam", -- uid of a medium defined under parameter or number
          saturation =
          { type = "vanGenuchten",
            value = "@SiltLoam",
          },
          conductivity =
          { type  = "vanGenuchten",
            value   = "@SiltLoam",
          },
      },
  },

  initial =
  {
    { cmp = "c", value = 0.0 },
    { cmp = "p", value = "HydroPressure" },
  },

  boundary =
  {
    -- Sea
    { cmp = "c", type = "dirichlet", bnd = "Sea", value = 1.0 },
    { cmp = "p", type = "dirichlet", bnd = "Sea", value = "HydroPressure_c" },

    -- Land
    { cmp = "c", type = "dirichlet", bnd = "Inflow", value = 0.0 },
    { cmp = "p", type = "flux", bnd = "Inflow", inner = "Medium", value = -3.3e-5*henry2D_rho}
  },

  linSolver =
  { type = "bicgstab",			-- linear solver type ["bicgstab", "cg", "linear"]
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
    dt		= 100,		            -- [s]  initial time step
    dtmin	= 0.00001,	          -- [s]  minimal time step
    dtmax	= tstop/10,	          -- [s]  maximal time step
    dtred	= 0.5,			          -- [1]  reduction factor for time step
    tol 	= 1e-2,
  },

  output =
  {
    file = "./", -- must be a folder!
    data = {"c", "p", "rho", "mu", "kr", "s", "q", "ff", "tf", "af", "df", "pc", "k"},
  },
}

function HydroPressure_c(x, y)
  return henry2D_g * 1025 * y
end

function HydroPressure(x, y)
  return henry2D_g * henry2D_rho * y
end

return henry
