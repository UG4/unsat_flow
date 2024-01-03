-- config for modelling a drainage trench with constant groundwater flow

local Qstrength = 2e+0

function TimeSink(x, y, t, si) 
  if (t<135.0) then return -Qstrength
  else return 0.0 end
end

-- command line parameters
params =
{
--	physical parameters
	fs_depth = util.GetParamNumber("-fsDepth", 0.25), -- depth of the free surface at the right boundary
	fs_slope = util.GetParamNumber("-fsSlope", 0), -- initial slope of the free surface
	
	recharge = util.GetParamNumber("-recharge", 0.0000165), -- "rain"
}

params.baseLvl = ARGS.numPreRefs

-- additional constants for vanGenuchten
rhog = 9.81 * 1000 

-- additional constants for vanGenuchten
local henry = 
{ 
  -- The domain specific setup
  domain = 
  {
    dim = 2,
    grid = "grids/henry_quad_sink.ugx",
    numRefs = ARGS.numRefs,
    numPreRefs = ARGS.numPreRefs,
  },

  -- list of non-linear models => translated to functions
  parameter = {  -- TODO: Parameters from List & Radu (2016)?
    { uid = "@Silt",
      type = "vanGenuchten",
      thetaS = 0.396, thetaR = 0.131,
      alpha = 0.423/rhog*10, n = 2.06, 
       }, 
    
    { uid = "@Clay",  -- modified n
      type = "vanGenuchten",
      alpha = 0.152/rhog*10, n = 3.06,  
      thetaS = 0.446, thetaR = 0.1, 
      },  --KSat= kappa/mu*rho*g   <=> kappa = Ksat*mu/(rho*g) 
    },

  flow = 
  {
    type = "haline",
    cmp = {"p", "c"},
    boussinesq = false,

    gravity = -9.81,    -- [ m s^{-2}�] ("standard", "no" or numeric value) 
    density =           
    { type = "linear",    -- density function ["linear", "exp", "ideal"]
      min = 1000, -- [ kg m^{-3} ] water density
      max = 1025.0,       -- [ kg m^{-3} ] saltwater density
      w_max = 1.0,
    },  
    
    viscosity = 
    { type = "const",      -- viscosity function ["const", "real"] 
      mu0 = 1e-3        -- [ kg m^{-3} ]  
    },
    diffusion   = 18.8571e-6,   -- constant
  },
   medium = 
   {
      {   subsets = {"Medium"}, 
          porosity = 0.396,
          saturation = 
          { type = "vanGenuchten",
            value = "@Silt",
          },
          conductivity =
          { type  = "vanGenuchten",
            value   = "@Silt",
          },
          permeability  = 1.019368e-9,  -- constant
      },
  },

  sources =
  {
     -- { cmp = "p", subset = "Sink", coord = {0.5, -0.75}, strength = -Qstrength},
     { 
        cmp = "p", subset = "Sink", coord = {0.5, -0.75}, strength = LuaUserNumber2d("TimeSink"), 
        substances = { 
          { cmp = "c" }, 
        },
     },
  },

  --[[sinks =
  {
     { 
        cmp = "c", subset = "Sink", coord = {0.5, -0.75}, strength = LuaUserNumber2d("TimeSink"),
        substances = {
         
        }
     },
  },
  ––]]


  initial = 
  {
    { cmp = "c", value = 0.0 },
    { cmp = "p", value = "HydroPressure" }
  },

  boundary = 
  {
    -- Sea
    { cmp = "c", type = "dirichlet", bnd = "Sea", value = 1.0 },
    { cmp = "p", type = "dirichlet", bnd = "Sea", value = "HydroPressure" },
    
    -- Land
    { cmp = "c", type = "dirichlet", bnd = "Inflow", value = 0.0 },
    { cmp = "p", type = "flux", bnd = "Inflow", inner = "Medium", value = -Qstrength },
    
    -- Recharge.
    -- { cmp = "p", type = "flux", bnd = "Top", inner="Medium", value="RechargeTop"}

  },

  solver =
  {
      type = "newton",
      lineSearch = {			   		-- ["standard", "none"]
          type = "standard",
          maxSteps		= 10,		-- maximum number of line search steps
          lambdaStart		= 1,		-- start value for scaling parameter
          lambdaReduce	= 0.5,		-- reduction factor for scaling parameter
          acceptBest 		= true,		-- check for best solution if true
          checkAll		= false		-- check all maxSteps steps if true 
      },

      convCheck = {
          type		= "standard",
          iterations	= 10,			-- number of iterations
          absolute	= 1e-8,			-- absolut value of defact to be reached; usually 1e-6 - 1e-9
          reduction	= 1e-7,		-- reduction factor of defect to be reached; usually 1e-6 - 1e-7
          verbose		= true			-- print convergence rates if true
      },
      
      linSolver =
      {
          type = "bicgstab",			-- linear solver type ["bicgstab", "cg", "linear"]
          precond = 
          {	
              type 		= "gmg",	-- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
              smoother 	= {type = "ilu", overlap = true},	-- gmg-smoother ["ilu", "ilut", "jac", "gs", "sgs"]
              cycle		= "V",		-- gmg-cycle ["V", "F", "W"]
              preSmooth	= 3,		-- number presmoothing steps
              postSmooth 	= 3,		-- number postsmoothing steps
              rap			= true,		-- comutes RAP-product instead of assembling if true 
              baseLevel	= params.baseLvl, -- gmg - baselevel
              
          },
          convCheck = {
              type		= "standard",
              iterations	= 30,		-- number of iterations
              absolute	= 0.5e-11,	-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
              reduction	= 1e-9,		-- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
              verbose		= true,		-- print convergence rates if true
          }
      }
  },
   
  time = 
  {
      control	= "limex",
      start 	= 0.0,				-- [s]  start time point
      stop	= 1000.0,			-- [s]  end time point
      max_time_steps = 10000,		-- [1]	maximum number of time steps
      dt		= ARGS.dt,		-- [s]  initial time step
      dtmin	= 0.00001 * ARGS.dt,	-- [s]  minimal time step
      dtmax	= 10.0,	-- [s]  maximal time step
      dtred	= 0.1,				-- [1]  reduction factor for time step
      tol 	= 1e-2,
  },
  
  output = 
  {
    freq	= 1, 	-- prints every x timesteps
    binary 	= true,	-- format for vtk file	
    file = "simulations/henry2D",
    data = {"c", "p", "q", "s", "k", "rho", "mu"}
  },
}

function HydroPressure_bnd(x, y, t, si) 
  pp = HydroPressure(x, y)
  if pp < 0 then return false, 0
  else 
  return true, pp 
  end
end

function HydroPressure(x, y) 
  return -10055.25 * (y + params.fs_depth) -- 9.81*1200
end

function RechargeTop(x, y, t, si) 
  return -(2.0-x)*3.3e-2
end



return henry

