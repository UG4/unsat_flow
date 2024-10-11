PrintBuildConfiguration()
-- ug_load_script("/home/software/d3f/ug4/apps/d3f_app/d3f_util.lua")

-----------------------------------------
-- Parameters --
-----------------------------------------
local YearInSeconds 	= 3.1557e+7 -- seconds
file_name = os.date("FG_Soil")
grid_name = "grids/soil2D.ugx"
diffu			= 1.0e-9
d_long			= 0 --.1
d_trans			= 0 --d_long/10
upwindmethod	= "full" -- "partial"
cp_steps		= 500
cp_name			= "CP_FG_Soil"
--T_stop			= 5e+2 * (YearInSeconds/12)
--deltaT			= T_stop/10000 --1e+3 * (YearInSeconds/12)
T_stop			= 5e+1 * (YearInSeconds/12)
deltaT			= T_stop/10000*10 --1e+3 * (YearInSeconds/12)
deltaT_min		= 1e-8 * (YearInSeconds/12)
deltaT_max		= 1e+5 * (YearInSeconds/12)
dt_start		= 3.1557e-2
freqoutput		= 1
--flow_vel		= 3e-8 	--m/s, wird für RB umgerechnet in kg/s
--rock_dens		= 2500

refs 			= 1		--refinement steps
prefs 			= 0 	--pre-refinement steps
lnslvty			= "gmg"
lnslvsmth		= "ilu"
preSmth			= 3
postSmth		= 3
newton_conv_abs = 1e-9
newton_conv_red = 1e-7
linSol_conv_abs	= 1e-10
linSol_conv_rel	= 1e-8


-- Arne (p=rho*g*h) 
-- rho * g 
--   1 Pa = N/m^2 = kg/m^3 * m/s^2 * m 
-- => m^2 / (Pa*s) = 
local CONST = {
  rho0 = 1000, 
  mu0 = 1e-3,
  g = -9.81,          -- [m/s^2]
  rhog = 9.81 * 1000, -- [Pa/m]
  DAY = 24*3600, -- [s]
  YEAR = 24*3600*360, -- [s]
}

-- Transformation:
--  kf = K/mu * rho *g            [m/s]
--  K  = (mu / (rho*g)) * kf      [m^2]

local function PermeabiltyFromConductivity(kf)
  local value = (CONST.mu0 / CONST.rhog) * kf 
  print("Permeability: "..value)
  return value
end

-- Equivalent of 0.25 m/a in kg/(m*m*s) 
local rechargeRate = 0.25 -- m/a / m^2
rechargeFlux = -(rechargeRate * CONST.rho0) / CONST.YEAR

-------------- ---------------------------
-- Functions
-----------------------------------------
-- hydrodynamic pressure
--function HydroPressuredepth_South(x,y,t) return (-y)*9810 end
function MyTracerSource(_,z,t)
	if t<=3.15571e+8 and z>=8 then s0=3.39e-06 else s0=0 end
	return s0
end

function MyHydrostaticP(x, z, t, si) 
  if (z<=5.0) then return (5.0-z) * CONST.rhog end-- 5.0 meters below ground level
  return -5000.0 -- (5.0-z) * CONST.rhog*0.0 end
end

local soil2D = 
{ 
  -- The domain specific setup
	domain = 
	{
		dim = 2,
		grid = grid_name, -- ("Ordner/File.ugx")
		numRefs = refs,
		numPreRefs = prefs,
	},

  -- list of non-linear models => translated to functions
  parameter = {  -- TODO: Parameters from List & Radu (2016)?
    
    { uid = "@MyVanGenuchten",
      type = "vanGenuchten",
     
      alpha = 1.4e-3,  
      n = 1.5, -- 1.5
      thetaS = 1.0, thetaR = 0.2,

    },--KSat= kappa/mu*rho*g   <=> kappa = Ksat*mu/(rho*g) 
    
  },
  
  flow = 
  {
    type = "haline",
    cmp = {"p", "c"},

    gravity = CONST.g,    -- [ m s^{-2}�] ("standard", "no" or numeric value) 
    density =           
    { 
      type = "const",         -- density function ["linear", "exp", "ideal"]
      min = CONST.rho0,       -- [ kg m^{-3} ] water density
    },  
    
    viscosity = 
    { type = "const",      -- viscosity function ["const", "real"] 
      mu0 = CONST.mu0        -- [ kg m^{-3} ]  
    },

    diffusion   = diffu, -- [ m^2/s ]
    upwind = upwindmethod

    --[[
    permeability 	=  	1.01936799184506E-12,
		alphaL			= 	d_long,
		alphaT			= 	d_trans,
		upwind 			= 	upwindmethod,	-- no, partial, full 
		boussinesq		= 	true,		-- true (no/less salt), false (with salt)/Einheiten Ausgabe Problem
	  
    datatable = 
		{
			{	"subset", 		"permeability", "porosity"},
			{	"soil", 		1.01936799184506E-11,	0.25},  
			{	"unsaturated",	1.01936799184506E-12, 	0.2	},
			{	"saturated", 	1.01936799184506E-12, 	0.2	}
		},
		
    --]]	
  },

  medium = 
  {
    {   
      subsets = {"soil"}, 
      porosity = 0.25,
      permeability  = PermeabiltyFromConductivity(1e-4), --  m^2,  -- constant m^2

      saturation = { value = "@MyVanGenuchten" },
      conductivity = { value = "@MyVanGenuchten"  },   -- relative permeabiltiy 
    
      -- Salt (=tracer) transport.
      diffusion   = diffu,   -- constant 
      alphaL      = d_long, 
      alphaT      = d_trans,  -- dispersion
    },

    {   
      subsets = {"unsaturated"}, 
      porosity = 0.2,
      permeability  = PermeabiltyFromConductivity(1e-5), --  m^2,  -- constant m^2

      saturation = { value = "@MyVanGenuchten" },
      conductivity = { value = "@MyVanGenuchten"  },   -- relative permeabiltiy 
    
      -- Salt (=tracer) transport.
      diffusion   = diffu,   -- constant 
      alphaL      = d_long, 
      alphaT      = d_trans,  -- dispersion
    },

    {   
      subsets = {"saturated"}, 
      porosity = 0.2,
      permeability  = PermeabiltyFromConductivity(1e-5), --  m^2,  -- constant m^2

      saturation = { value = "@MyVanGenuchten" },
      conductivity = { value = "@MyVanGenuchten"  },   -- relative permeabiltiy 
    
      -- Salt (=tracer) transport.
      diffusion   = diffu,   -- constant 
      alphaL      = d_long, 
      alphaT      = d_trans,  -- dispersion
    },
  },

  sources = {},

  initial = 
   {
       { cmp = "p", value = "MyHydrostaticP" },
       { cmp = "c", value = "MyTracerSource" },
   },

  boundary = 
  {
     -- Fresh water influx on TOP
     {cmp = "p", type = "flux", bnd = "Top", inner="soil", value = rechargeFlux},
     {cmp = "c", type = "dirichlet", bnd = "Top", value = 0.0},
     
     -- pressure Dirichlet BC (equivalent to 5m water column)  
     {cmp = "p", type = "dirichlet", bnd = "Base", value = CONST.rhog * 5.0 },
     {cmp = "c", type = "dirichlet", bnd = "Base", value = 0.0},

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
          iterations	= 128,			-- number of iterations
          absolute	= 1e-8,			-- absolut value of defact to be reached; usually 1e-6 - 1e-9
          reduction	= 1e-7,		-- reduction factor of defect to be reached; usually 1e-6 - 1e-7
          verbose		= true			-- print convergence rates if true
      },
      
      linSolver =
      {
          type = "bicgstab",			-- linear solver type ["bicgstab", "cg", "linear"]
          precond = 
          {	
              type 		= lnslvty,	-- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
              smoother 	= {type = "ilu", overlap = true},	-- gmg-smoother ["ilu", "ilut", "jac", "gs", "sgs"]
              cycle		= "V",		-- gmg-cycle ["V", "F", "W"]
              preSmooth	= preSmth,		-- number presmoothing steps
              postSmooth 	= postSmth,		-- number postsmoothing steps
              rap			= true,		-- comutes RAP-product instead of assembling if true 
              baseLevel	= 0, -- gmg - baselevel
              
          },
          convCheck = {
              type		= "standard",
              iterations	= 100,		-- number of iterations
              absolute	= linSol_conv_abs,	-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
              reduction	= linSol_conv_rel,		-- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
              verbose		= true,		-- print convergence rates if true
          }
      }
  },
  time = 
  {
      control	= "limex",
      start 	= 0.0,				-- [s]  start time point
      stop	= 10.0*CONST.YEAR,			-- [s]  end time point
      max_time_steps = 1000,		-- [1]	maximum number of time steps
      dt		= 1e-4*ARGS.dt*CONST.DAY,		-- [s]  initial time step
      dtmin	= 1e-14 * ARGS.dt*CONST.DAY,	-- [s]  minimal time step
      dtmax	= 10.0*CONST.DAY,	-- [s]  maximal time step
      dtred	= 0.1,				-- [1]  reduction factor for time step
      tol 	= 1e-3,

      -- startstep 			= {dt = dt_start},
		  -- checkpoint_steps 	= cp_steps,
		  -- checkpoint_name 	= cp_name,
  },
  
  --[[
  TODO: Include this code

  flux =
  {
  {name = "Top",			bnd = "Top", 	file="Top.dat"},
  {name = "Base",			bnd = "Base", 	file="Base.dat"},
  {name = "Sides",		bnd = "Sides", 	file="Sides.dat"},		
  {name = "Top_Tracer",	bnd = "Top", 	file="Top_Tracer.dat", cmp = "Tracer"},
  {name = "Base_Tracer",	bnd = "Base", 	file="Base_Tracer.dat", cmp = "Tracer"},
  {name = "Sides_Tracer",	bnd = "Sides", 	file="Sides_Tracer.dat", cmp = "Tracer"},
  },

output = 
{
  vtkname 				= file_name,	-- name of vtk file
      write_pvd 				= true,
      write_processwise_pvd 	= false,
  binary 					= true,
  freq					= freqoutput,
  vtktimes 				= 
    {
      1.0*YearInSeconds/12,
      10.0*YearInSeconds/12,
      20.0*YearInSeconds/12,
      30.0*YearInSeconds/12,
      40.0*YearInSeconds/12,
      50.0*YearInSeconds/12,
      60.0*YearInSeconds/12,
      70.0*YearInSeconds/12,
      80.0*YearInSeconds/12,
      90.0*YearInSeconds/12,
      100.0*YearInSeconds/12,
      200.0*YearInSeconds/12,
      300.0*YearInSeconds/12,
      400.0*YearInSeconds/12,
      500.0*YearInSeconds/12,
    },
  {file = "Output_Top_BC_flux.dat",		type = "flux", 	data ="q", boundary="Top", 		inner= "soil"},
  {file = "Output_Base_BC_flux.dat",		type = "flux", 	data ="q", boundary="Base", 	inner= "saturated"},
  {file = "Output_Sides_BC_flux.dat", 	type = "flux", 	data ="q", boundary="Sides", 	inner= "soil,unsaturated,saturated"},
  {file = "Tracer_Integral_all.dat", 		type = "integral", data = "Tracer"},
}
--]]
  -- config for vtk output
  -- possible data output variables: 
  -- c (concentration), p (pressure), q (Darcy Velocity), s (saturation),
  -- k (conductivity), f (flux), rho (density)
  output = 
    {
      freq	= 1, 	-- prints every x timesteps
      binary 	= true,	-- format for vtk file	
      file = "soil2D",
      data = {"c", "p", "q", "s", "vf", "rho", "kr"},
      
      fs_evaluation_points = { 
        {0.0}
      }
    }
}




return soil2D

