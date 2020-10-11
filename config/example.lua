problem = 
{ 
	-- The domain specific setup
	domain = 
	{
		dim = 2,                        -- problem dimension
		grid = "grids/levee2D.ugx",    -- grid name
		numRefs = 2,                    -- refinements after distribution
		numPreRefs = 1,                 -- refinements before distribution
	},

	-- The density-driven-flow setup
	flow = 
	{
		type = "haline",
		cmp = {"c", "p"},

		gravity = -9.81,    -- [ m s^{-2}�] ("standard", "no" or numeric value)	
		density = 					
		{	type = "linear", 		-- density function ["linear", "exp", "ideal"]
			min = 1000,				-- [ kg m^{-3} ] water density
			max = 1200				-- [ kg m^{-3} ] saltwater density
		},	
		
		viscosity = 
		{	type = "real",			-- viscosity function ["const", "real"] 
			mu0 = 1e-3				-- [ kg m^{-3} ]	
		},
		air_pressure = 1013.25e2,
	},

	medium = 
	{
	    { 	subsets = {"CLAY"}, 
	      	porosity = 1.0,
			saturation = 
			{	type = "vanGenuchten",
				value = "@Silt"
			},
	      	conductivity =
			{	type	= "vanGenuchten", 
				value 	= "@Silt"
			},
			diffusion		= 3.565e-6, 	-- constant
			permeability 	= 4.845e-13, 	-- constant

	    },
	    { 	subsets = {"SAND_LEFT","SAND_RIGHT"}, 
	      	porosity = 1.0,
			saturation 		= 
			{	type = "vanGenuchten",
				value = "@Silt"
			},
	      	conductivity	=
			{	type	= "vanGenuchten",
				value	= "@Silt"
			},
			diffusion		= 3.565e-6,
			permeability 	= 4.845e-13,
	    },
	},

	parameter = 
	{
	    { uid = "@Silt", 
	      type = "vanGenuchten",
	      thetaS = 0.396, thetaR = 0.131,
	    --thetaS = 0.396, thetaR = 0.396,
	      alpha = 0.423, n = 2.06, Ksat= 4.96e-1 }, 
	    
	    { uid = "@Clay",  -- modified n
	      type = "vanGenuchten",
	      alpha = 0.152, n = 3.06,  
	      thetaS = 0.446, thetaR = 0.1, 
	      Ksat= 8.2e-4 * 1e-3,}, 
    },

	initial_conditions = 
	{
	    { cmp = "p", value = "ConcentrationStart" },
	},

	boundary_conditions = 
	{
	   	{cmp = "p", type = "dirichlet", bnd = "WaterBnd", value = "ConcentrationDirichletBnd"},
	    {cmp = "p", type = "dirichlet", bnd = "ToeBnd", value = 0.0 },
	},
		
	solver =
	{
		type = "newton",
		lineSearch = {			   		-- ["standard", "none"]
			type = "standard",
			maxSteps		= 30,		-- maximum number of line search steps
			lambdaStart		= 1,		-- start value for scaling parameter
			lambdaReduce	= 0.5,		-- reduction factor for scaling parameter
			acceptBest 		= true,		-- check for best solution if true
			checkAll		= false		-- check all maxSteps steps if true 
		},

		convCheck = {
			type		= "standard",
			iterations	= 100,			-- number of iterations
			absolute	= 5e-8,			-- absolut value of defact to be reached; usually 1e-6 - 1e-9
			reduction	= 1e-20,		-- reduction factor of defect to be reached; usually 1e-6 - 1e-7
			verbose		= true			-- print convergence rates if true
		},
		
		linSolver =
		{
			type = "bicgstab",			-- linear solver type ["bicgstab", "cg", "linear"]
			precond = 
			{	
				type 		= "gmg",	-- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
				smoother 	= "ilu",	-- gmg-smoother ["ilu", "ilut", "jac", "gs", "sgs"]
				cycle		= "V",		-- gmg-cycle ["V", "F", "W"]
				preSmooth	= 3,		-- number presmoothing steps
				postSmooth 	= 3,		-- number postsmoothing steps
				rap			= false,	-- comutes RAP-product instead of assembling if true 
				baseLevel	= 0,		-- gmg - baselevel
				
			},
			convCheck = {
				type		= "standard",
				iterations	= 60,		-- number of iterations
				absolute	= 1e-10,	-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
				reduction	= 1e-3,		-- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
				verbose		= true,		-- print convergence rates if true
			}
		}
	},
	
	time = 
	{
		control	= "prescribed",
		start 	= 0.0,			-- [s]  start time point
		stop	= 3.1536e8,		-- [s]  end time point
		dt		= 3.1536e6,		-- [s]  initial time step
		dtmin	= 3.1536e3,		-- [s]  minimal time step
		dtmax	= 3.1536e6,		-- [s]  maximal time step
		dtred	= 0.1,			-- [1]  reduction factor for time step
		tol 	= 1e-2,
		
	},
	
	output = 
	{
		freq	= 1, 			-- prints every x timesteps
		binary 	= true,			-- format for vtk file
		vtkname = "Elder",		-- name of vtk file
		
	}
} 

-- Some functions (sadly global...)

Levee2D_tRise = levee2D.var.RiseTime

-- Anstieg bis auf 5.85 m
function Levee2D_Pegel(t) 
  return math.min(t/Levee2D_tRise, 1.0)*5.85 
end  

-- This is a rising flood.
function Levee2D_RisingFlood(x, z, t, si) 
  local zPegel = Levee2D_Pegel(t)
  if (z <zPegel) then return true, math.max(zPegel-z, 0.0) end
  return false, 0.0
end

-- This is an empty levee.
function Levee2D_HydrostaticHead(x, z, t, si) 
  return -z
end

return levee2D




--function ConcentrationStart(x, y, t)
--	if 