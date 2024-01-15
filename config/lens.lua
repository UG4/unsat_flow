-- config for modelling a drainage trench with constant groundwater flow

well_rho = 998.23
well_g = -9.81 -- must be negative!
rhog = (-1.0)*well_rho*well_g
numdays = 1000
tstop = numdays * 86400 -- 100 days

local lens =
{
  -- The domain specific setup
  domain =
  {
    dim = 2,
    grid = "grids/lens.ugx",
    numRefs = ARGS.numRefs,
    numPreRefs = ARGS.numPreRefs,
  },

  -- medium parameters for vanGenuchten Model
  parameter = {
    { uid = "@Sandstone",
      type = "vanGenuchten",
      thetaS = 0.250, thetaR = 0.153,
      alpha = 0.79/rhog, n = 10.4,
      Ksat = 1.08/86400},

    { uid = "@TouchetSiltLoam",
      type = "vanGenuchten",
      thetaS = 0.469, thetaR = 0.190,
      alpha = 0.50/rhog, n = 7.09,
      Ksat = 3.03/86400},

    { uid = "@SiltLoam",
      type = "vanGenuchten",
      thetaS = 0.396, thetaR = 0.131,
      alpha = 0.423/rhog, n = 2.06,
      Ksat = 0.0496/86400},
  },

  flow =
  {
    boussinesq = false,

    gravity = well_g,      -- [ m s^{-2}], must be negative!
    density =
    { type = "ideal",         -- density function ["linear", "exp", "ideal"]
      min = well_rho,      -- [ kg m^{-3} ] water density
      max = 1025.0,           -- [ kg m^{-3} ] saltwater density
    },
    diffusion   = 18.8571e-6 -- [ m^2/s ]
  },
  medium =
  {
     {   subsets = {"Inner"},
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
    {cmp = "p", strength = "pump", subset = "Well", coord = {1.0, 1.5}, substances = {cmp = "c"}}
  },

  initial =
  {
    { cmp = "c", value = 1.0 },
    { cmp = "p", value = "WellPressureStart" },
  },

  boundary =
  {
    -- Top
    {cmp = "p", type = "flux", bnd = "Top", inner="Inner", value = -0.001},
    {cmp = "c", type = "dirichlet", bnd = "Top", value = 0.0},

    -- Sea
    {cmp = "p", type = "dirichlet", bnd = "Sea", value = "WellSeaBoundary" },
    {cmp = "c", type = "dirichlet", bnd = "Sea", value = 1.0},

    -- Bottom
    {cmp = "p", type = "dirichlet", bnd = "Bottom", value = "WellSeaBoundary"},
    {cmp = "c", type = "dirichlet", bnd = "Bottom", value = 1.0},
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
    dt		= 1,		          -- [s]  initial time step
    dtmin	= 0.001,	          -- [s]  minimal time step
    dtmax	= tstop/10,	            -- [s]  maximal time step
    dtred	= 0.5,			          -- [1]  reduction factor for time step
    tol 	= 1e-2,
  },

  output =
  {
    file = "./", -- must be a folder!
    data = {"c", "p", "rho", "mu", "kr", "s", "q", "ff", "tf", "af", "df", "pc", "k"},
  }

}


function WellSeaBoundary(x, y, t)
  return true, (1.0 - y) * rhog
end

function WellPressureStart(x, y, t)
  return (1.0 - y) * rhog
end

function pump(x, y, t)
  if (t > 10*86400) then
    return -0.001
  else
    return 0.0
  end
end

return lens
