-- Modelled after:
-- Improving the worthiness of the Henry problem as a benchmark for density-dependent groundwater flow models
-- by Matthew J. Simpson and T. Prabhakar Clement
-- 2004

-- command line parameters
params =
{
--	physical parameters
	fs_depth = util.GetParamNumber("-fsDepth", 0.25), -- depth of the free surface at the right boundary
	recharge = util.GetParamNumber("-recharge", 0.0000165), -- "rain"
  frac_rain = util.GetParamNumber("-frac", 0.0), -- "rain"
}

-- additional constants for vanGenuchten
henry2D_rho = 1000
henry2D_g = -9.8 -- must be negative!
rhog = (-1.0)*henry2D_rho*henry2D_g
-- tstop = 1000 -- seconds -- 100 * 86400 -- 100 days
tstop = 5000 -- 100 * 86400 -- 100 days


local henry =
{
  -- The domain specific setup
  domain =
  {
    dim = 2,
    -- grid = "grids/henry_quad.ugx",
    grid = "grids/henry_tri_4.ugx",
    numRefs = ARGS.numRefs,
    numPreRefs = ARGS.numPreRefs,
  },

  -- medium parameters for vanGenuchten Model
  parameter = {
    { uid = "@Sand",
      type = "vanGenuchten",
      thetaS = 1.0, thetaR = 0.0,
      alpha = 0.087/rhog*1e+0,
      n = 1.58,
      Ksat = 1e-2} -- [ m/s ]
  },

  flow =
  {
    boussinesq = false,
    
    gravity = henry2D_g,      -- [ m s^{-2}], must be negative!
    density =
    { type = "linear",         -- density function ["linear", "exp", "ideal"]
      min = henry2D_rho,      -- [ kg m^{-3} ] water density
      max = 1025.0,           -- [ kg m^{-3} ] saltwater density
    },
    diffusion   = 18.8571e-6,  -- [ m^2/s ]
   -- diffusion   = 0.0,  -- [ m^2/s ]
    upwind  = "full"
  },
   medium =
   {
      {   subsets = {"Medium"},
          porosity = 0.35, -- uid of a medium defined under parameter or number
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
    { cmp = "c", type = "dirichlet", bnd = "Sea", value = "SeaConc" },
    { cmp = "p", type = "dirichlet", bnd = "Sea", value = "HydroPressure_bnd" },

    -- Land -- Closed
    --{ cmp = "c", type = "dirichlet", bnd = "Inflow", value = 0.0 },
    --{ cmp = "p", type = "flux", bnd = "Inflow", inner = "Medium", value = -3.3e-2*}

    -- Top
    { cmp = "p", type = "flux", bnd = "Top", inner="Medium", value=-3.3e-2*0.5*0.0},
    { cmp = "c", type = "dirichlet", bnd = "Top", value = 0.0 },

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
     --[[ convCheck =
      { type    = "composite",
        iterations  = 30,   -- number of iterations
        absolute  = 0.5e-8, -- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
        reduction = 1e-7,   -- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
        verbose   = true ,   -- print convergence rates if true
        sub ={
        {cmp ="p"},
        {cmp="c"}
        }
      }]]--
  },

  time =
  {
    control	= "limex",
    start 	= 0.0,				      -- [s]  start time point
    stop	= tstop,			        -- [s]  end time point
    max_time_steps = 1000,		  -- [1]	maximum number of time steps
    dt		= tstop*1e-6,		            -- [s]  initial time step
    dtmin	= tstop*1e-9,	        -- [s]  minimal time step
    dtmax	= tstop/100,	            -- [s]  maximal time step
    dtred	= 0.5,			          -- [1]  reduction factor for time step
    tol 	= 1e-3,
  },

  output =
  {
    file = "./", -- must be a folder!
    data = {"c", "p", "rho", "mu", "kr", "s", "q", "ff", "tf", "af", "df", "pc", "k", "vf"},
  },
}

function HydroPressure_c(x, y)
  return henry2D_g * 1025 * (y + params.fs_depth)
end

function HydroPressure_bnd(x, y, t, si)
  local pp = HydroPressure_c(x, y)
  if pp < 0 then return false, 0 end
  return true, pp
end

function SeaConc(x, y, t, si)
  local pp = HydroPressure_c(x, y)
  if pp < 0 then return false, 0.0 end
  return true, 1.0 
end

function HydroPressure(x, y)
  return henry2D_g * henry2D_rho * (y + params.fs_depth)
end

function RechargeTop(x, y, t, si)
  return -(2.0-x)*params.recharge
end


return henry
