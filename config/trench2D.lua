-- config for modelling a drainage trench with constant groundwater flow

Trench2D_rho = 998.23
Trench2D_g = -9.81 -- must be negative!
rhog = (-1.0)*Trench2D_rho*Trench2D_g

local trench2D =
{
  -- The domain specific setup
  domain =
  {
    dim = 2,
    grid = "grids/trench2D.ugx",
    numRefs = ARGS.numRefs,
    numPreRefs = ARGS.numPreRefs,
  },

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

    gravity = Trench2D_g, -- [m s^{-2}]
    density =
    { type = "ideal",     -- density function ["linear", "exp", "ideal"]
      min = Trench2D_rho, -- [ kg m^{-3} ] water density
      max = 1025.0,       -- [ kg m^{-3} ] saltwater density
    },

    viscosity =
    { type = "real",      -- viscosity function ["const", "real"]
      mu0 = 1.002e-3      -- [ kg m^{-3} ]
    },
  },
   medium =
   {
      {   subsets = {"Inner"},
          porosity = "@SiltLoam",
          saturation =
          { type = "vanGenuchten",
            value = "@SiltLoam",
          },
          conductivity =
          { type  = "vanGenuchten",
            value   = "@SiltLoam",
          },
          diffusion   = 18.8571e-6,   -- constant
          permeability  = "@SiltLoam" -- 1.019368e-9,  -- must be uid of a medium defined under parameter or number
      },
  },

  initial=
   {
       { cmp = "p", value = "Trench2DPressureStart"},
       { cmp = "c", value = 0}
   },

  boundary =
  {
     {cmp = "p", type = "dirichlet", bnd = "Trench", value = "Trench2DDrainagePressureBoundary"},
     {cmp = "p", type = "dirichlet", bnd = "Aquifer", value = "Trench2DAquiferBoundary" },
     {cmp = "c", type = "dirichlet", bnd = "Trench", value = 1},
     {cmp = "c", type = "dirichlet", bnd = "Aquifer", value = 0},
  },

  solver =
  {
      type = "newton",
      lineSearch = {			   		  -- ["standard", "none"]
          type = "standard",
          maxSteps		= 10,		    -- maximum number of line search steps
          lambdaStart		= 1,		  -- start value for scaling parameter
          lambdaReduce	= 0.5,    -- reduction factor for scaling parameter
          acceptBest 		= true,   -- check for best solution if true
          checkAll		= false		  -- check all maxSteps steps if true
      },

      convCheck = {
          type		= "standard",
          iterations	= 128,    -- number of iterations
          absolute	= 1e-8,	    -- absolut value of defact to be reached; usually 1e-6 - 1e-9
          reduction	= 1e-7,		  -- reduction factor of defect to be reached; usually 1e-6 - 1e-7
          verbose		= true	    -- print convergence rates if true
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
              iterations	= 30,		-- number of iterations
              absolute	= 0.5e-8,	-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
              reduction	= 1e-7,		-- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
              verbose		= true,		-- print convergence rates if true
          }
      }
  },

  time =
  {
    control = "limex",
    start   = 0.0,          -- [s] start time point
    stop  = 2000.0,         -- [s] end time point
    dt  = 0.01,             -- [s] initial time step
    max_time_steps = 10000, -- [1]	maximum number of time steps
    dtmin	= 0.01 * ARGS.dt,	-- [s]  minimal time step
    dtmax	= 10.0,	          -- [s]  maximal time step
    dtred = 0.1,            -- [1] reduction factor for time step
    tol   = 1e-2
  },

  output =
  {
    file = "simulations/trench2D/", -- needs to be a folder!
    data = {"c", "p", "rho", "mu", "kr", "s", "q", "f", "pc", "k"}
  }

}


function Trench2DDrainagePressureBoundaryTime(x, y, t, tD)
  if (t <= tD) then
    return true, (2.2*t / tD - 2.0) * (-1.0) * rhog
  else
    return true, 0.2  * (-1.0) * rhog
  end
end

function Trench2DDrainagePressureBoundary(x, y, t)
  return Trench2DDrainagePressureBoundaryTime(x, y, t, 220.0)
end

function Trench2DAquiferBoundary(x, y, t)
  return true, (1.0 - y) * rhog
end

function Trench2DPressureStart(x, y, t)
  return (1.0 - y) * rhog
end

return trench2D
