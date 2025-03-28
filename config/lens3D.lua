local rho0 = 1000 -- kg/m^3
local rhog = 9.81 * 1000 -- approx: 1e+4 -- kg/(m*s)
local rhog1 = 9.81 *1020
local z0 = 5.0 -- meters below ground level.
local alpha = 1e+1
local Ss = 0.25 -- 1/m


local DAY = 3600*24 -- seconds
local myRadius = 100.0 -- 

local myHeight = 10.0

local myRechargeRate = 0.1/DAY -- m/day -- Recharge Rate
local myRechargeRadius = 12.5 -- 





local mySectorArea = (4.0/3.0*math.pi*myRadius*myRadius)*0.125
local myRechargeArea = (4.0/3.0*math.pi*myRechargeRadius*myRechargeRadius)*0.125  -- m^2
local myRechargeVol = myRechargeArea*myRechargeRate -- cbm/s
local myRechargeMass = myRechargeVol*rho0 -- kg/s


local QStrength = 0.5 * myRechargeMass

local lens3D = 
{ 
  -- The domain specific setup
  domain = 
  {
    dim = 3,
    grid = "grids/lens3D_d20.ugx",
    numRefs = ARGS.numRefs,
    numPreRefs = ARGS.numPreRefs,
    -- neededSubsets = {}
  },

  -- list of non-linear models => translated to functions
  parameter = { 
    
    { uid = "@Silt",
      type = "vanGenuchten",
      thetaS = 0.396, thetaR = 0.131,
      alpha = 0.423/rhog, n = 2.06},--4.96e-1 -- }, 
    
    { uid = "@Clay",  -- modified n
      type = "vanGenuchten",
      alpha = 0.152/rhog, n = 3.06,  
      thetaS = 0.446, thetaR = 0.1},  --KSat= kappa/mu*rho*g   <=> kappa = Ksat*mu/(rho*g) 

    { uid = "@MyExponential", 
      type = "exp",
      pentry = rhog,
      alpha = 1e+1, -- TODO: Storativity?
      beta = 0.0, -- 0.0, --2e+1, 
      thetaR = 1e-6, thetaS = 0.1}, 
    
        --[[
    { uid = "@WaltherSaturation", 
      type = "const",
      value = 1.0
    }, 
    --]]
    -- [[
    { 
    	-- S(p)  = thetaS * exp(alpha*p/pentry)
    	-- S'(p) = thetaS * alpha/pentry S(p) \approx thetaS * alpha/pentry
    	uid = "@WaltherSaturation", 
      type = "exp",
      pentry = rhog, -- normalize to height
      alpha = alpha, 			-- TODO: Storativity?
      beta = 0.0, -- 0.0, --2e+1, 
      thetaR = 1e-12, thetaS = 1.0/alpha 
    }, --]]
    
   --[[
    { uid = "@WaltherPermeabilty", 
      type = "const",
      value = 1.0 -- Relative permeability
    }, 
    --]]
    

  },
  
  flow = 
  {
    type = "haline",
    cmp = {"p", "c"},

    gravity = -9.81,    -- [ m s^{-2}�] ("standard", "no" or numeric value) 
    density =           
    { type = "linear",    -- density function ["linear", "exp", "ideal"]
      min = 1000,       -- [ kg m^{-3} ] water density
      max = 1020,       -- [ kg m^{-3} ] saltwater density
      w_max = 1,
    },  
    
    viscosity = 
    { type = "const",      -- viscosity function ["const", "real"] 
      mu0 = 1e-3        -- [ kg m^{-3} ]  
    },

    diffusion   = 1e-9, -- [ m^2/s ]
    upwind = "partial"
  },

  medium = 
  {
    {   
      subsets = {"Aquifer"}, 
      porosity = 0.2,
      
      saturation = { value = "@Silt" },
      conductivity = { value = "@Silt"  }, -- relative permeabiltiy
      

      --conductivity = { value = "@WaltherSaturation" }, -- relative permeabiltiy
      --saturation = {  value = "@WaltherSaturation" },
      --saturation = { ["$ref"] = "./well3D.lua#/parameters/WaltherSaturation" } 
      --saturation = { ["$ref"] = "#/parameters/WaltherSaturation" } 
      -- cf. https://redocly.com/docs/resources/ref-guide
     
      --saturation = {  value = 1.0 },
      --conductivity = { value = 1.0 }, -- relative permeabiltiy
      -- storativity = { value = 5.0968E-08 }, 
      --mass_storage = { value = rho0*Ss/rhog}, -- /9.81}, -- relative permeabiltiy
    
      permeability  = 1.0194e-12, --  m^2 (from Walther)
        
      -- Salt transport.
      diffusion   = 1.0e-3,   -- constant 
      alphaL = 0.0, alphaT = 0.0,           -- w/o dispersion
    },
  },

  sources =
  {
    
     -- [[ 
     { 
        cmp = "p", subset = "Sink", coord = {0.0, 0.0, -20.0}, 
        strength = -QStrength, 
        substances = 
        { 
          { cmp = "c" },
        },
     },
     --]]
  },

  initial = 
   {
       { cmp = "p", value = "Lens_3D_Initial_P" },
       { cmp = "c", value = "Lens_3D_Initial_C"},
   },

  boundary = 
  {
     -- {cmp = "p", type = "flux", inner = "Aquifer", bnd = "Rim", value = -QStrength/mySectorArea},
    {cmp = "p", type = "dirichlet", bnd = "Rim", value = "Lens_3D_Seaside_BC_P"},
    {cmp = "c", type = "dirichlet", bnd = "Rim", value = "Lens_3D_Seaside_BC_C"},
     
     
     {cmp = "p", type = "flux", inner = "Aquifer", bnd = "Top", value = "Lens_3D_Top_BC_P"},
     {cmp = "c", type = "dirichlet", bnd = "Top", value = "Lens_3D_Top_BC_C"},
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
              type 		= "gmg",	-- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
              smoother 	= {type = "ilu", overlap = true},	-- gmg-smoother ["ilu", "ilut", "jac", "gs", "sgs"]
              cycle		= "V",		-- gmg-cycle ["V", "F", "W"]
              preSmooth	= 3,		-- number presmoothing steps
              postSmooth 	= 3,		-- number postsmoothing steps
              rap			= true,		-- comutes RAP-product instead of assembling if true 
              baseLevel	= 0, -- gmg - baselevel
              baseSolver=SuperLUCPU1()
              
          },
     
          convCheck = {
            type		= "composite",
            iterations	= 100,		-- number of iterations
            absolute	= 0.5e-8,	-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
            reduction	= 1e-7,		-- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
            verbose		= true,		-- print convergence rates if true
            sub = 
            {
              {cmp = "c", absolute = 1e-12, relative=1e-8},
              {cmp = "p", absolute = 1e-12, relative=1e-9}
            }
        }
        --]]  
        
      }
  },
  time = 
  {
      control	= "limex",
      start 	= 0.0,				-- [s]  start time point
      stop	= 20*360.0*DAY,			-- [s]  end time point
      max_time_steps = 100000,		-- [1]	maximum number of time steps
      dt		= 1e-0*ARGS.dt*DAY,		-- [s]  initial time step
      dtmin	= 1e-14 * ARGS.dt*DAY,	-- [s]  minimal time step
      dtmax	= 360.0*DAY,	-- [s]  maximal time step
      dtred	= 0.1,				-- [1]  reduction factor for time step
      tol 	= 1e-3,

      metricSpace = {
        H1SemiComponentSpace3dCPU1("c", 2),
        L2ComponentSpace3dCPU1("p", 2)
      },

      
      -- Idea: new items (JSON reference):
      post_process = { ["$ref"] = "./lens3D.lua#/TheisToolbox/PostProcess" } -- JSON-like reference to LUA function.
  },
  
  
  -- config for vtk output
  -- possible data output variables: 
  -- c (concentration), p (pressure), q (Darcy Velocity), s (saturation),
  -- k (conductivity), f (flux), rho (density)
  output = 
  {
    freq	= 1, 	-- prints every x timesteps
    binary 	= true,	-- format for vtk file	
    file = "simulations/lens3D",
    data = {"c", "p", "q", "s", "kr", "rho"},
      
    fs_evaluation_points = { 
    -- [[
      {0.0, 0.0}, {3.125, 0.0}, {9.375, 0.0}, {12.0, 0.0},
      {18.75, 0.0}, {25.0, 0.0}, {37.5, 0.0}, {50.0, 0.0},
      {75.0, 0.0}, {100.0, 0.0}
    --]]
    } 
  },
  
  -- 
  functions =
  {
    -- We can put arbitrary functions here. 
    HydrostaticPressure = function (x, y, z, t, si) 
      return (-z0-z) * rhog1 -- z0 meters below ground level
    end

  }
}


-- Some functions (unfortunately as global variables...)
function Lens_3D_Initial_P(x, y, z, t, si) 
  return (-z0-z) * rhog1 -- z0 meters below ground level (saltwater)
end

function Lens_3D_Initial_C(x, y, z, t, si) 
  if z<= -z0 then
    return 1.0 -- z0 meters below ground level
  end
   return 0.0
end

function Lens_3D_Seaside_BC_P(x, y, z, t, si) 
  --if z<= -z0 then
   return true, (-z0-z) * rhog1 -- z0 meters below ground level
  --end
  --return true, 0
end

function Lens_3D_Seaside_BC_C(x, y, z, t, si) 
  if z<= -z0 then return true, 1.0 -- z0 meters below ground level
  end
  return true, 0.0
end


function Lens_3D_Top_BC_P(x, y, z, t, si) 
  local r2=x*x+y*y
  if (r2<=myRechargeRadius*myRechargeRadius) then return -1.0*myRechargeRate*rho0  -- recharge in m/s
  else return 0
  end
end

function Lens_3D_Top_BC_C(x, y, z, t, si) 
  local r2=x*x+y*y
  --if (r2<=myRechargeRadius*myRechargeRadius) then 
    return true, 0.0 -- z0 meters below ground level
  --else return false, 0.0
  --end   
end
 

print("==============================================================")
-- print("V="..0.125*mySectorVol)
print("Q="..QStrength)

local myMedium = lens3D.medium[1]

local myStorativity = Ss*(10-z0)
print("Storativity (S)    [1]      =\t\t"..myStorativity)

local myConductivity = (myMedium.permeability/ lens3D.flow.viscosity.mu0)* rhog
local myTransmissivity =  myConductivity*(10-z0)
print("Transmissivity (T) [m*m/s]  =\t\t" .. myTransmissivity)
print("Conductivity (Kf)  [m/s]    =\t\t" .. myConductivity)
print("Recharge Rate (R)  [m/s]    =\t\t" .. myRechargeRate)
print("Ratio (R/Kf)       [1]      =\t\t" .. myRechargeRate/myConductivity)



--myStorativity = myMedium.mass_storage.value * rhog
--print("Storativity (S2)   =\t\t"..myStorativity)
print("S/T                 [s/m*m] =\t\t"..(myStorativity/myTransmissivity))
print("S/T  (Walther)       [s/m*m] =\t\t25.000")
print("Ss/Kf               [s/m*m] =\t\t"..(Ss/myConductivity))
print("Dhyd=Kf/Ss          [m*m/s]=\t\t"..(myConductivity/Ss))

-- print("ss/(kappa/mu)      =\t\t"..(Ss/myConductivity))

print("Char. time diff     [s] =\t\t"..(100*100*Ss/myConductivity))
print("Char. time pump     [s] =\t\t".. 1.0/(QStrength/rho0))

TheisToolbox = {
  -- Store value
  ST=(Ss/myConductivity),
  Q = QStrength*8.0/rho0,  -- analytic source 
  T = myTransmissivity,
}
-- Theis well function approximation in Lua
function theis_well_function(u, nterms)
  -- Check if the input is valid
  if u <= 0 then
      --error("Parameter 'u' must be greater than 0")
      return 0.0/0.0
  end

  nterms = nterms or 10 -- Default number of terms in the series expansion
  local gamma = 0.57721566490153286060 -- Euler-Mascheroni constant
  local result = -math.log(u) - gamma
  local term = u
  local factorial = 1

  --print(result)
  -- Series expansion loop

  if (u<10.0) then
    for n = 1, nterms do
      local update = term / (factorial*n)
      result = result + update
      term = term * u * (-1.0)
      factorial = factorial * (n + 1)
      --print(update)
      if math.abs(update) / math.abs(result) < 1e-16 then
         return result 
      end
    end
  return result
  else 
     -- Use asymptotic expansion for large u
     local result = 0
     local term = 1 / u
     local factorial = 1
     local exp_neg_u = math.exp(-u)

     for n = 0, nterms - 1 do
         result = result + term
         term = -term * (n + 1) / u
     end

     return exp_neg_u * result
    end
end


TheisToolbox.compute_ref = function (x,y,t,si)
    local u = (x*x+y*y)/(4.0*t) * TheisToolbox.ST
    local w = theis_well_function(u, 20)
    return w * TheisToolbox.Q / (4.0* math.pi * TheisToolbox.T)
end

-- Evaluate Theis' well solution (for a given set of points).
TheisToolbox.post_process = function (points, step, time, currdt)

  -- local points = lens3D.output.fs_evaluation_points
  -- print(points) -- DEBUG: print evaluation points
  print("-------- Theis (t="..time..")---------")
  local dd ={}
 
   for i, x in ipairs(points) do
    local res=TheisToolbox.compute_ref(x[1], x[2], time, 0)
    --print(res)
    dd[i] = res

    -- print("{"..x[1].. "," ..x[2].."} \t".. time .. "\t".. res)
    print(string.format("{%.3e, %.3e} \t %.4e \t %e " , x[1], x[2], time, -z0-res))
  end
  
  -- Print drawdown
  -- print(dd)
end

-- print(TheisToolbox)

-- Define PostProcess (TODO: functions break JSON compatibility...)
lens3D.PostProcess = function (step, time, currdt) 

  -- Evaluate reference
  TheisToolbox.post_process(lens3D.output.fs_evaluation_points, step, time, currdt)


end

print("==============================================================")


return lens3D

