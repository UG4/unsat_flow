-- Freshwater lens benchmark modelled after Stoeckl et al. A new numerical benchmark of a freshwater lens (2016)

lens_rho = 997
lens_rho_c = 1021
lens_g = -9.81 -- must be negative!
rhog = (-1.0)*lens_rho*lens_g
tstop = 24 * 60 * 60 -- 1 day

local lens =
{
  -- The domain specific setup
  domain =
  {
    dim = 2,
    grid = "grids/stoeckl_lens.ugx",
    numRefs = ARGS.numRefs,
    numPreRefs = ARGS.numPreRefs,
  },

  -- medium parameters for vanGenuchten Model
  parameter = {
    { uid = "@Material",
      type = "vanGenuchten",
      thetaS = 0.39, thetaR = 0.1,
      alpha = 0.423/rhog, n = 2.06,
      Ksat = 4.1e-3}
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
    upwind = "full"
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
         },
         alphaL = 5e-4,
         alphaT = 0.1*5e-4,
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
    {cmp = "p", type = "flux", bnd = "Top", inner="Inner", value = "top_boundary"},
    {cmp = "c", type = "dirichlet", bnd = "Top", value = "top_boundary_c"},

    -- Left
    {cmp = "p", type = "dirichlet", bnd = "Left", value = "left_boundary"},
    {cmp = "c", type = "dirichlet", bnd = "Left", value = 1.0},
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
    dtmax	= 8.64*10,	            -- [s]  maximal time step
    dtred	= 0.3,			          -- [1]  reduction factor for time step
    tol 	= 1e-2,
  },

  output =
  {
    file = "./", -- must be a folder!
    data = {"c", "p", "rho", "mu", "kr", "s", "q", "ff", "tf", "af", "df", "pc", "k"},
  }

}

T0 = 6*60*60 -- phase switch, after quasi steady state, 6h

function HydroPressure(x, y)
  return y * lens_rho_c * lens_g
end

function top_boundary(x, y, t, si)
  -- Shoreline segment of 1cm between 0.5m and 0.51m
  if t >= T0 or x <= 0.51 then
    return false, 0.0
  else
    return true, -1.333e-5*lens_rho -- [m^3/s]
  end
end

function top_boundary_c(x, y, t, si)
  if t >= T0 or x <= 0.51 then
    return false, 0.0
  else
    return true, 0.0
  end
end


function left_boundary(x, y, t, si)
  return true, y * lens_rho_c * lens_g 
end

return lens
