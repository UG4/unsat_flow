-- Freshwater lens benchmark modelled after Stoeckl et al. A new numerical benchmark of a freshwater lens (2016)

lens_rho = 997
lens_rho_c = 1021
lens_g = -9.81 -- must be negative!
rhog = (-1.0)*lens_rho*lens_g

recharge_rate = util.GetParamNumber("--recharge", -1.8e-5)
total_time = util.GetParamNumber("--hours", 24, "Total simulation time in hours")
steady_state = 60*60*10 -- 10 hours to reach steady state
pump_rate = -5e-3 -- 1.0 m^3/day
sea_level = util.GetParamNumber("--sea_level", 0.27, "Sea level in m") -- 0.3 for fully saturated

tstop = total_time * 60 * 60

function HydroPressure(x, y)
  return (y - sea_level) * lens_rho_c * lens_g -- phreatic surface at y = 0.4m
end

function pumping(x, y, t, si)
  if t < steady_state then
    return 0.0
  else
    return pump_rate
  end
end

function shore_boundary(x, y, t, si)
  hp = HydroPressure(x, y)
  if y > sea_level then
    return false, 0.0 -- no flow above sea level
  else
    return true, hp
  end
end

function shore_boundary_c(x, y, t, si)
  if y > sea_level then
    return false, 0.0
  else
    return true, 1.0
  end
end

local lens =
{
  -- The domain specific setup
  domain =
  {
    dim = 2,
    grid = "grids/stoeckl_lens_pump_full.ugx",
    numRefs = ARGS.numRefs,
    numPreRefs = ARGS.numPreRefs,
  },

  -- medium parameters for vanGenuchten Model
  parameter = {
    { uid = "@Material",
      type = "vanGenuchten",
      thetaS = 0.39, thetaR = 0.1,
      alpha = 0.423/rhog, n = 2.06,
      Ksat = 4.5e-3}
  },

  flow =
  {
    boussinesq = false,

    gravity = lens_g,      -- [ m s^{-2}], must be negative!
    density =
    { type = "linear",         -- density function ["linear", "exp", "ideal"]
      min = lens_rho,      -- [ kg m^{-3} ] water density
      max = lens_rho_c,           -- [ kg m^{-3} ] saltwater density
    },
    diffusion   = 10e-9, -- [ m^2/s ]
    upwind = "partial"
  },
  medium =
  {
     {   subsets = {"Inner"},
         porosity = "@Material", -- uid of a medium defined under parameter or number
         saturation =
         { type = "vanGenuchten",
           value = "@Material",
         },
         conductivity =
         { type  = "vanGenuchten",
           value   = "@Material",
         }
     },
 },

 initial =
  {
    { cmp = "c", value = 1.0 },
    { cmp = "p", value = "HydroPressure" },
  },

  boundary =
  {
    -- Top
    {cmp = "p", type = "neumann", bnd = "Top", inner="Inner", value = recharge_rate*lens_rho },
    {cmp = "c", type = "dirichlet", bnd = "Top", value = 0.0 },

    -- Shoreline
    {cmp = "p", type = "dirichlet", bnd = "Shore", value = "shore_boundary"},
    {cmp = "c", type = "dirichlet", bnd = "Shore", value = "shore_boundary_c"},
  },

  sources =
  {
    {cmp = "p", strength = "pumping", subset = "Pump", coord = {0.9, 0.15}, substances = {{cmp = "c"}}},
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
    max_time_steps = 10000,		  -- [1]	maximum number of time steps
    dt		= 0.864,		          -- [s]  initial time step
    dtmin	= 0.001,	          -- [s]  minimal time step
    dtmax	= 8.64*1000,	            -- [s]  maximal time step
    dtred	= 0.3,			          -- [1]  reduction factor for time step
    tol 	= 1e-2,
  },

  output =
  {
    file = "./", -- must be a folder!
    data = {"c", "p", "rho", "mu", "kr", "s", "q", "ff", "tf", "af", "df", "pc", "k"},
    fs_evaluation_points = { 
       {0.0, 0.0}, {0.1, 0.0}, {0.2, 0.0}, {0.3, 0.0}, {0.4, 0.0}, {0.5, 0.0}, {0.6, 0.0}, {0.7, 0.0}, {0.8, 0.0}, {0.9, 0.0}
     } 
  }

}

return lens
