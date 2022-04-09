------------------------------------------------------------------------------------------------
-- Lua script for the simulation of the Stoeckl setting.
--
-- Cf. L. Stoeckl, M. Walther, L. K. Morgan, Physical and Numerical Modelling of Post-Pumping Seawater Intrusion, DOI: 10.1155/2019/7191370
------------------------------------------------------------------------------------------------


tstop = 300 * 86400 -- 100 days
params =
{
}

params.freshWaterDensity = 996.9
params.saltWaterDensity = 1020.9
params.gravity = -9.981

rhog = (-1.0)*params.freshWaterDensity*params.gravity

local henry =
{
  -- The domain specific setup
  domain =
  {
    dim = 2,
    grid = "grids/ext_stoeckl_quad.ugx",
    numRefs = ARGS.numRefs,
    numPreRefs = ARGS.numPreRefs,
  },

  -- medium parameters for vanGenuchten Model
  parameter = {
    { uid = "@SiltLoam",
      type = "vanGenuchten",
      thetaS = 0.41, thetaR = 0.131,
      alpha = 0.423, n = 2.06,
      Ksat = 2.2e-3}
  },

  flow =
  {
    boussinesq = false,

    gravity = params.gravity,      -- [ m s^{-2}], must be negative!
    density =
    { type = "linear",         -- density function ["linear", "exp", "ideal"]
      min = params.freshWaterDensity,	-- [ kg m^{-3} ]
      max = params.saltWaterDensity,	-- [ kg m^{-3} ]
    },

    viscosity =
    { type = "const",          -- viscosity function ["const", "real"]
      mu0 = 0.00089                  -- [ Pa s ]
    },
    diffusion   = 10e-9, -- [ m^2/s ]
    upwind = "partial"
  },
   medium =
   {
      {   subsets = {"inner"},
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

  sources =
  {
    {cmp = "p", value = "Pumping", subset = "pump", x = 0.6, y = 0.05},
    {cmp = "c", transport = "Pumping", subset = "pump", x = 0.6, y = 0.05},
  },

  initial =
  {
    { cmp = "c", value = 0.0 },
    { cmp = "p", value = "HydroPressureInitial" },
  },

  boundary =
  {
    -- Sea
    { cmp = "c", type = "dirichlet", bnd = "salt", value = 1.0 },
    { cmp = "p", type = "dirichlet", bnd = "salt", value = "HydroPressure" },

    -- Land
    { cmp = "p", type = "flux", bnd = "fresh", inner = "inner", value = -1.11e-5 }, -- approx. 0.96 [m3/d]
    { cmp = "c", type = "dirichlet", bnd = "fresh", value = 0.0 }
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
    dt		= 0.2,		          -- [s]  initial time step
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

function HydroPressure(x,y)
  return params.gravity * params.saltWaterDensity * (y - 0.12) end

function HydroPressureInitial(x,y)
  return params.gravity * params.saltWaterDensity * (y - 0.2) end


function Pumping(x, y, t, si)
  t_phase1 = 5000 -- first phase, 5000s
  t_phase2 = 135  -- second phase, 135s
  t_phase3 = 5000 -- third phase, 5000s
  if t < t_phase1 then
    return 0
  elseif t > t_phase1 and t < (t_phase1 + t_phase2 + t_phase3) then
    return -1.157e-5 -- 1 m/d
  end
  return 0
end

return henry
