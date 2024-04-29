------------------------------------------------------------------------------------------------
-- Lua script for the simulation of the Stoeckl setting.
--
-- Cf. L. Stoeckl, M. Walther, L. K. Morgan, Physical and Numerical Modelling of Post-Pumping Seawater Intrusion, DOI: 10.1155/2019/7191370
------------------------------------------------------------------------------------------------

params =
{
}

params.freshWaterDensity = 996.9
params.saltWaterDensity = 1020.9
params.gravity = -9.981

params.pumping_rate = -1 / 86400 * params.freshWaterDensity -- 1 [m3/d] spread along the entire width
params.phase1_time = 5000 -- 5000s: steady state generation
params.phase2_time = 135-- 135s: pumping
params.phase3_time = 5000 -- 5000s: pump shut off

tstop = 35000 -- 12 hours

-- The fresh water recharge through the right boundary:
-- 0.086 [m3/d] (cf. the paper, p. 3, line 18 from below), converted to [kg/s] and
-- spread along the entire height (so that the integral provides the same mass)
params.freshflux = -0.96 / 86400 * params.freshWaterDensity -- [kg/s]

function HydroPressure(x, y, t, si)
	return params.gravity * params.saltWaterDensity * y
end

function pumping_sink(x, y, t, si)
	if t < params.phase1_time then
    return 0.0
  elseif t < params.phase2_time then
    return params.pumping_rate
  else
    return 0.0
  end
end

local henry =
{
  -- The domain specific setup
  domain =
  {
    dim = 2,
    grid = "grids/stoeckl_pump.ugx",
    numRefs = ARGS.numRefs,
    numPreRefs = ARGS.numPreRefs,
  },

  -- medium parameters for vanGenuchten Model
  parameter = {
    { uid = "@SiltLoam",
      type = "vanGenuchten",
      thetaS = 0.41, thetaR = 0.131,
      alpha = 0.423, n = 2.06,
      Ksat = 1.5e-3}
  },

  flow =
  {
    gravity = params.gravity,      -- [ m s^{-2}], must be negative!
    density =
    { type = "linear",         -- density function ["linear", "exp", "ideal"]
      min = params.freshWaterDensity,	-- [ kg m^{-3} ]
      max = params.saltWaterDensity,	-- [ kg m^{-3} ]
    },
    diffusion   = 1e-9, -- [ m^2/s ]
    upwind = "full",
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
          alphaL = 5e-3,
          alphaT = 5e-4,
      },
  },

  --[[sources =
  {
    {cmp = "p", strength = "pumping_sink", subset = "pump", coord = {0.6, 0.05}, substances = {{cmp = "c"}}},
  },--]]

  initial =
  {
    { cmp = "c", value = 0.0 },
    { cmp = "p", value = "HydroPressure" },
  },

  boundary =
  {
    -- Left
    { cmp = "c", type = "dirichlet", bnd = "Salt", value = 1.0 },
    { cmp = "p", type = "dirichlet", bnd = "Salt", value = "HydroPressure" },

    -- Right
    { cmp = "p", type = "flux", bnd = "Inflow", inner = "Medium", value = params.freshflux },
    { cmp = "c", type = "dirichlet", bnd = "Inflow", value = 0.0 }
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
    dtmin	= 0.0001,	          -- [s]  minimal time step
    dtmax	= 2500,	          -- [s]  maximal time step
    dtred	= 0.5,			          -- [1]  reduction factor for time step
    tol 	= 1e-2,
  },

  output =
  {
    file = "./", -- must be a folder!
    data = {"c", "p", "rho", "mu", "kr", "s", "q", "ff", "tf", "af", "df", "pc", "k"},
  },
}


return henry
