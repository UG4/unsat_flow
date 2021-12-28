-- config for modelling the henry problem

-- command line parameters
params =
{
--	physical parameters
	fs_depth = util.GetParamNumber("-fsDepth", 0.0), -- depth of the free surface at the right boundary
	recharge = util.GetParamNumber("-recharge", 0.0000165), -- "rain"
}

-- additional constants for vanGenuchten
henry2D_rho = 998.23
henry2D_g = -9.81 -- must be negative!
rhog = (-1.0)*henry2D_rho*henry2D_g
tstop = 50 * 86400 -- 100 days


local henry =
{
  -- The domain specific setup
  domain =
  {
    dim = 2,
    grid = "grids/henry_quad.ugx",
    numRefs = ARGS.numRefs,
    numPreRefs = ARGS.numPreRefs,
  },

  -- medium parameters for vanGenuchten Model
  parameter = {
    { uid = "@Sand",
      type = "vanGenuchten",
      thetaS = 0.37, thetaR = 0.043,
      alpha = 0.087/rhog, n = 1.58,
      Ksat = 1}
  },

  flow =
  {
    boussinesq = false,

    gravity = henry2D_g,      -- [ m s^{-2}], must be negative!
    density =
    { type = "ideal",         -- density function ["linear", "exp", "ideal"]
      min = henry2D_rho,      -- [ kg m^{-3} ] water density
      max = 1025.0,           -- [ kg m^{-3} ] saltwater density
    },

    viscosity =
    { type = "const",          -- viscosity function ["const", "real"]
      mu0 = 1.002e-3                   -- [ Pa s ]
    },
    diffusion   = 7.64e-8 -- [ m^2/s ]
  },
   medium =
   {
      {   subsets = {"Medium"},
          porosity = "@Sand", -- uid of a medium defined under parameter or number
          saturation =
          { type = "vanGenuchten",
            value = "@Sand",
          },
          conductivity =
          { type  = "vanGenuchten",
            value   = "@Sand",
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
    { cmp = "p", type = "dirichlet", bnd = "Sea", value = "HydroPressure_bnd" },

    -- Land
    { cmp = "c", type = "dirichlet", bnd = "Inflow", value = 0.0 },
    { cmp = "p", type = "flux", bnd = "Inflow", inner = "Medium", value = -7.64e-5 }

    -- Top
    --{ cmp = "p", type = "flux", bnd = "Top", inner="Medium", value=params.recharge},
    --{ cmp = "c", type = "dirichlet", bnd = "Top", value = 0.0 },

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
    --[[convCheck =
      { type		= "standard",
        iterations	= 30,		-- number of iterations
        absolute	= 0.5e-8,	-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
        reduction	= 1e-7,		-- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
        verbose		= true		-- print convergence rates if true
      }]]--
      convCheck =
      { type    = "composite",
        iterations  = 30,   -- number of iterations
        absolute  = 0.5e-8, -- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
        reduction = 1e-7,   -- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
        verbose   = true ,   -- print convergence rates if true
        sub ={
        {cmp ="p"},
        {cmp="c"}
        }
      } 
  },

  time =
  {
    control	= "limex",
    start 	= 0.0,				      -- [s]  start time point
    stop	= tstop,			        -- [s]  end time point
    max_time_steps = 1000,		  -- [1]	maximum number of time steps
    dt		= 1000,		          -- [s]  initial time step
    dtmin	= ARGS.dt,	          -- [s]  minimal time step
    dtmax	= tstop/10,	            -- [s]  maximal time step
    dtred	= 0.5,			          -- [1]  reduction factor for time step
    tol 	= 1e-2,
  },

  output =
  {
    file = "./", -- must be a folder!
    data = {"c", "p", "rho", "mu", "kr", "s", "q", "ff", "tf", "af", "df", "pc", "k"},
  },
}

function HydroPressure_bnd(x, y, t, si)
  pp = HydroPressure_c(x, y)
  if pp < 0 then
    return false, 0
  else
    return true, pp
  end
end

function HydroPressure_c(x, y)
  return henry2D_g * 1025 * (y + params.fs_depth)
end

function HydroPressure(x, y)
  return henry2D_g * henry2D_rho * (y + params.fs_depth)
end

function RechargeTop(x, y, t, si)
  return -(2.0-x)*params.recharge
end


return henry
