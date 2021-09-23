-- config for modelling the henry problem

-- command line parameters
params =
{
--	physical parameters
	fs_depth = util.GetParamNumber("-fsDepth", 0.2), -- depth of the free surface at the right boundary
	recharge = util.GetParamNumber("-recharge", 3.3e-2), -- "rain"
}

-- additional constants for vanGenuchten
henry2D_rho = 998.23
henry2D_g = -9.81 -- must be negative!
rhog = (-1.0)*henry2D_rho*henry2D_g

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
    { uid = "@Sandstone",
      type = "vanGenuchten",
      thetaS = 0.153, thetaR = 0.250,
      alpha = 0.79/rhog, n = 10.4,
      Ksat = 1.08},

    { uid = "@TouchetSiltLoam",
      type = "vanGenuchten",
      thetaS = 0.190, thetaR = 0.469,
      alpha = 0.50/rhog, n = 7.09,
      Ksat = 3.03},

    { uid = "@SiltLoam",
      type = "vanGenuchten",
      thetaS = 0.396, thetaR = 0.131,
      alpha = 0.423/rhog, n = 2.06,
      Ksat = 0.0496},

    { uid = "@SiltLoam",
      type = "vanGenuchten",
      thetaS = 0.446, thetaR = 0.0,
      alpha = 0.152/rhog, n = 1.17,
      Ksat = 8.2e-4}
    },

  flow =
  {
    boussinesq = false,

    gravity = henry2D_g,      -- [ m s^{-2}], must be negative!
    density =
    { type = "ideal",     -- density function ["linear", "exp", "ideal"]
      min = henry2D_rho,         -- [ kg m^{-3} ] water density
      max = 1025.0,       -- [ kg m^{-3} ] saltwater density
    },

    viscosity =
    { type = "real",      -- viscosity function ["const", "real"]
      mu0 = 1.002e-3         -- [ kg m^{-3} ]
    },
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
          diffusion   = 18.8571e-6,   -- constant
          permeability  = = "@SiltLoam" -- 1.019368e-9,  -- must be uid of a medium defined under parameter or number
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
    -- { cmp = "p", type = "flux", bnd = "Inflow", inner = "Medium", value = -3.3e-2 },
    { cmp = "p", type = "flux", bnd = "Top", inner="Medium", value="RechargeTop"}

  },

  solver =
  {
      type = "newton",
      lineSearch = {			   		  -- ["standard", "none"]
          type = "standard",
          maxSteps		= 10,		    -- maximum number of line search steps
          lambdaStart		= 1,		  -- start value for scaling parameter
          lambdaReduce	= 0.5,		-- reduction factor for scaling parameter
          acceptBest 		= true,		-- check for best solution if true
          checkAll		= false		  -- check all maxSteps steps if true
      },

      convCheck = {
          type		= "standard",
          iterations	= 10,			-- number of iterations
          absolute	= 1e-8,			-- absolut value of defact to be reached; usually 1e-6 - 1e-9
          reduction	= 1e-7,		  -- reduction factor of defect to be reached; usually 1e-6 - 1e-7
          verbose		= true			-- print convergence rates if true
      },

      linSolver =
      {
          type = "bicgstab",			-- linear solver type ["bicgstab", "cg", "linear"]
          precond =
          {
              type 		= "gmg",	                          -- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
              smoother 	= {type = "ilu", overlap = true},	-- gmg-smoother ["ilu", "ilut", "jac", "gs", "sgs"]
              cycle		= "V",		                          -- gmg-cycle ["V", "F", "W"]
              preSmooth	= 3,		                          -- number presmoothing steps
              postSmooth 	= 3,		                        -- number postsmoothing steps
              rap			= true,		                          -- comutes RAP-product instead of assembling if true
              baseLevel	= ARGS.numPreRefs,                -- gmg - baselevel

          },
          convCheck = {
              type		= "standard",
              iterations	= 30,		  -- number of iterations
              absolute	= 0.5e-11,	-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
              reduction	= 1e-9,		  -- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
              verbose		= true,		  -- print convergence rates if true
          }
      }
  },

  time =
  {
      control	= "limex",
      start 	= 0.0,				      -- [s]  start time point
      stop	= 1000.0,			        -- [s]  end time point
      max_time_steps = 10000,		  -- [1]	maximum number of time steps
      dt		= ARGS.dt,		        -- [s]  initial time step
      dtmin	= 0.00001 * ARGS.dt,	-- [s]  minimal time step
      dtmax	= 10.0,	              -- [s]  maximal time step
      dtred	= 0.1,			          -- [1]  reduction factor for time step
      tol 	= 1e-2,
  },

  output =
  {
    file = "simulations/henry2D", -- needs to be a folder!
    data = {"c", "p", "rho", "mu", "kr", "s", "q", "f", "pc", "k"}
  },
}

function HydroPressure_bnd(x, y, t, si)
  pp = HydroPressure(x, y)
  if pp < 0 then
    return false, 0
  else
    return true, pp
  end
end

function HydroPressure(x, y)
  return -10055.25 * (y + params.fs_depth)
end

function RechargeTop(x, y, t, si)
  return -(2.0-x)*recharge
end


return henry
