-- Lua-Script to perform an unsaturated density flow problem
-- Author: Niklas Conen

local myPath = ug_get_current_path()
package.path = package.path .. ";" .. myPath .. "/config/?.lua;" .. myPath .. "/?.lua"
-- package.path = package.path .. ";" .. myPath .. "/config/trench/?.lua"
package.path = package.path .. ";./config/?.lua;./?.lua"
print(package.path)

ug_load_script("../scripts/ug_util.lua")
ug_load_script("./unsat_flow_util.lua")
ug_load_script("../scripts/util/solver_util.lua")

util.CheckAndPrintHelp("Unsaturated density flow problem");

-- Parameters every problem uses
-- problem specific parameters are in the config file
ARGS =
{
  problemID  = util.GetParam("--problem-id", "trench2D"),
  numPreRefs = util.GetParamNumber("--numPreRefs", 1, "number of refinements before parallel distribution"),
  numRefs    = util.GetParamNumber("--numRefs", 4, "number of refinements after parallel distribution"),
  adaptive   = util.HasParamOption("--adaptive", false),
  dt         = util.GetParamNumber("--dt", 0.01, "time step"),
  validate   = util.HasParamOption("--validate", true),
}

print (ARGS.problemID)
function unsatSolve(problemID, numPreRefs, numRefs, adaptive)
  if adaptive then
    ug_load_script("./unsat_flow_adaptive.lua")
  end

  local problem

  -- if the problemID is given as string, load the problem from the config file
  -- if it is given as object, use the object
  if type(problemID) == "string" then
    problem = require(problemID)
  else
    problem = problemID
  end

  vtools = require("validation")
  if (validate and vtools and vtools.validate(problem)) then
    print ("Problem '".. problemID .. "' validated successfully!")
  else
    print ("Problem '".. problemID .. "' is invalid (or cannot be validated)!")
    -- exit() --remove this line, if you prefer to continue with 'invalid' schemes...
    print(problem)
  end

  InitUG(problem.domain.dim, AlgebraType("CPU", 1))

  local dom = util.CreateAndDistributeDomain(problem.domain.grid, numRefs, numPreRefs, {})

  -- saves the refined grid
  -- SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), "refined.ugx", 0.1)

  local disc = ProblemDisc:new(problem, dom)

  -- create approximation space.
  local approxSpace = disc:CreateApproxSpace()

  disc.u = GridFunction(disc.approxSpace)

  -- Creating the Domain discretisation for the problem
  local domainDisc = disc:CreateDomainDisc(approxSpace)

  -- vtk output
  disc.vtk = VTKOutput()
  disc:CreateVTKOutput()
  print("Created VTK Output")

  -- Initial Data
  disc:SetInitialData(disc.u)

  -- Solver Config
  util.solver.defaults.approxSpace = approxSpace

  -- Time stepping parameters
  local startTime                  = problem.time.start
  local endTime                    = problem.time.stop
  local dt                         = problem.time.dt
  local dtMin                      = problem.time.dtmin
  local dtMax                      = problem.time.dtmax
  local TOL                        = problem.time.tol
  local dtred                      = problem.time.dtred

  local dbgWriter                  = GridFunctionDebugWriter(approxSpace)

  -- LIMEX time-stepping
  -- Solvers config.
  local nstages = 2
  local limexLSolver = {}
  local limexNLSolver = {}
  local limexDomainDisc = {}

  if problem.linSolver == nil then
    problem.linSolver = problem.solver.linSolver
  end

  -- Linear solvers.
  for i = 1, nstages do
    limexLSolver[i] = util.solver.CreateSolver(problem.linSolver)
    limexDomainDisc[i] = domainDisc
  end

  if (not adaptive) then
    -- One-step Newton
    local limexConvCheck = ConvCheck(1, 1e-12, 1e-10, true)
    limexConvCheck:set_supress_unsuccessful(true)
    limexConvCheck:set_minimum_defect(3e-10)

    for i = 1, nstages do
      limexNLSolver[i] = NewtonSolver()
      limexNLSolver[i]:set_linear_solver(limexLSolver[i])
      limexNLSolver[i]:set_convergence_check(limexConvCheck)
    end
  else
    local niConfig = {
      maxLvl = 8, -- largest level
      maxSteps = 3,
      rTOL = 1e-4,  -- spatial tolerance
      aTOL = 1e-12, -- spatial tolerance

      refiner = HangingNodeDomainRefiner(dom),

      spaceOPT = H1SemiComponentSpace("p", 2, 1.0),
      markerOPT = GradientIndicator(disc.u, "p"),

      --refinementMarking = MaximumMarking(0.8, 0.0, 0.03),   -- max, min, fraction discarded
      --refinementMarking = EquilibrationMarkingStrategy(0.9,0.1),
      refinementMarking = VarianceMarkingEta(1, 0.8, 0.05),
      --refinementMarking = VarianceMarkingEta(3, 0.8, 0.05),
      coarseningMarking = APosterioriCoarsening(0.025), --VarianceMarkingEta(0.25, 3, 0.25),-- MaximumMarking(1.0, 0.05, 0.0),

      -- Debug.
      useVTK = true, -- true for debug
    }

    adaptivityUtil = AdaptivityUtil:new(niConfig, approxSpace, domainDisc, domainDisc)
    print(adaptivityUtil)

    -- Nested iteration solvers.
    local niSolverAdaptive = adaptivityUtil:CreateNestedIterationSolver(limexLSolver[1])
    limexNLSolver[1] = niSolverAdaptive;
    limexDomainDisc[1] = domainDisc

    -- Error for debug.
    local approxSpace_pc = ApproximationSpace2d(dom, AlgebraType("CPU", 1))
    approxSpace_pc:add_fct("eta_squared", "piecewise-constant");
    -- approxSpace_pc:add_fct("dummy", "piecewise-constant");


    local dbgEta = GridFunction2dCPU1(approxSpace_pc)
    dbgEta:set(0.0)
    niSolverAdaptive:set_debug_elem_error(dbgEta)
    niSolverAdaptive:set_debug(dbgWriter)


    for i = 2, nstages do
      local niSolverNonAdaptive = adaptivityUtil:CreateNestedIterationSolver(limexLSolver[i])
      niSolverNonAdaptive:disable_adaptive_refinement()
      limexNLSolver[i] = niSolverNonAdaptive
      limexDomainDisc[i] = domainDisc
    end

    print(limexDomainDisc)
  end

  -- LIMEX Configuration
  --- descriptor for LIMEX time integrator
  local limexDesc = {
    nstages = nstages,
    steps = { 1, 2, 3, 4, 5, 6 },

    domainDisc = limexDomainDisc,
    nonlinSolver = limexNLSolver,

    tol = 3e-3,
    dt = problem.time.dt,
    dtmin = problem.time.dtmin,
    dtmax = problem.time.dtmax,

    rhoSafetyOPT = 0.8,
    spacesOPT = WeightedH1SemiSpace,
    --debugOPT = dbgWriter,
  }

  print("Creating LIMEX integrator...")
  print(limexDesc)
  local limex = util.limex.CreateIntegrator(limexDesc)
  print("...DONE!")

  local weightedMetricSpace = CompositeSpace()

  if false then
    local scaleP = 1.0
    local spaceP = H1SemiComponentSpace("p", 2)
 
    local scaleC = 1e+12
    local spaceC = L2ComponentSpace("c", 2, scaleC)

    weightedMetricSpace:add(spaceC, scaleC)
    weightedMetricSpace:add(spaceP, scaleP)
  else
    -- Scale  with || (kappa_0/mu_0) * grad(p) || 
    -- Scale  with || (kappa_0/mu_0) * rho' * g  * w || 

    local kappa_over_mu_squared = 1.0 -- 4.60095884e-7 
    local spaceP = VelEnergyComponentSpace("p", 2, ConstUserMatrix(kappa_over_mu_squared))
    local spaceC = L2ComponentSpace("c", 2, kappa_over_mu_squared*(200*10)*200*10)

    weightedMetricSpace:add(spaceP)
    weightedMetricSpace:add(spaceC)
  end

  local limexErrorEst = CompositeGridFunctionEstimator()
  limexErrorEst:add(weightedMetricSpace)

  limex:add_error_estimator(limexErrorEst)
  limex:set_tolerance(problem.time.tol)
  limex:set_stepsize_safety_factor(0.8)
  limex:set_time_step(problem.time.dt)
  limex:set_dt_min(problem.time.dtmin)
  limex:set_dt_max(problem.time.dtmax)
  limex:set_increase_factor(10.0)
  --limex:enable_matrix_cache()
  limex:disable_matrix_cache() -- recompute solution.

  limex:set_conservative(true)

  -- Debugging LIMEX.

  if (false) then
    limex:set_debug(dbgWriter)
    limex:set_debug_for_timestepper(dbgWriter)
    limexNLSolver[1]:set_debug(dbgWriter)
  end

  if adaptive then
    print(adaptivityUtil)
    local luaobserver = LuaCallbackObserver()

    -- work-around (waiting for implementation of SmartPtr forward to lua...)
    function luaAdaptivePostProcess(step, time, currdt)
      print("This is step " .. step)
      adaptivityUtil:GridCoarseningStep(luaobserver:get_current_solution(),
        step, time, currdt)
      return 0;
    end

    limex:attach_observer(luaobserver)
    luaobserver:set_callback("luaAdaptivePostProcess")
  end


  -- Debugging LIMEX.
  local dbgWriter = GridFunctionDebugWriter(approxSpace)
  if (false) then
    --limex:set_debug(dbgWriter)
    --limex:set_debug_for_timestepper(dbgWriter)
    dbgWriter:set_vtk_output(true)
    dbgWriter:set_conn_viewer_output(true)
    limexNLSolver[1]:set_debug(dbgWriter)
  end

  -- Time step observer.
  local filename = (problem.output.filename) or (problemID)
  if (filename == nil) then
    filename = "unsat_flow"
  end

  local vtkobserver = nil
  if problem.output.plot_step then
    vtkobserver = VTKOutputObserver(problem.output.file..filename..".vtk", disc.vtk, problem.output.plot_step)
  else
    vtkobserver = VTKOutputObserver(problem.output.file..filename..".vtk", disc.vtk)
  end

  limex:attach_observer(vtkobserver)

  -- post process for saving time step size
  local luaObserver = LuaCallbackObserver()
  function luaPostProcessIntegration(step, time, currdt)
    print("")
    print("INT_V00 - Integral over salt mass fraction:\t" .. time .. "\t" .. Integral(disc.CompositeVoluFrac, disc.u))
    print("INT_V01 - Integral over salt mass fraction:\t" .. time .. "\t" .. Integral(disc.CompositeSaltMass, disc.u))
    print("INT_V02 - Integral over fluid phase volume:\t".. time .. "\t" .. Integral(disc.CompositeSaturation, disc.u))
    --print("INT_S02a - Integral gradp*n (in):\t" .. time .. "\t" .. IntegrateNormalGradientOnManifold(disc.u, "p","Inflow", "Medium" ))
    --print("INT_S02b - Integral gradp*n (out):\t" .. time .. "\t".. IntegrateNormalGradientOnManifold(disc.u, "p","Sea", "Medium" ))

    --print("INT_S03a - Integral gradc*n (in):\t" .. time .. "\t" .. IntegrateNormalGradientOnManifold(disc.u, "c","Inflow", "Medium" ))
    --print("INT_S03b - Integral gradc*n (out):\t" .. time .. "\t".. IntegrateNormalGradientOnManifold(disc.u, "c","Sea", "Medium" ))
    print("")
    print(">>>> TimeStep: " .. step .. "," .. time .. "," .. currdt .. " <<<<")
    print("")
    return 0;
  end

  luaObserver:set_callback("luaPostProcessIntegration")
  limex:attach_observer(luaObserver)

  ----------------------------------------
  -- BEGIN phreatic surface observer.
  ----------------------------------------
  local myFSPoints = problem.output.fs_evaluation_points
  if (myFSPoints and type(myFSPoints) == "table" and #myFSPoints>0
      and FSFileMeasurer) then
    surfaceObserver = FSFileMeasurer(disc.u, "p", 0.0)
    surfaceObserver:enable_print_output()
    surfaceObserver:enable_step_file_output("phreatic_surface_height")
    saltInterfaceObserver = FSFileMeasurer(disc.u, "c", 0.5)
    saltInterfaceObserver:enable_print_output()
    saltInterfaceObserver:enable_step_file_output("salt_interface_height")

    -- Add points.
    for k,vpoint in pairs(myFSPoints) do
      surfaceObserver:add_measurement_point(vpoint)
      saltInterfaceObserver:add_measurement_point(vpoint)
    end

    limex:attach_observer(surfaceObserver)
    limex:attach_observer(saltInterfaceObserver)
  end

  ----------------------------------------
  -- END phreatic surface observer.
  ----------------------------------------

  local sw = CuckooClock()
  sw:tic()
  -- Solve problem.

  limex:apply(disc.u, endTime, disc.u, 0.0)

  print("CDELTA=" .. sw:toc())
end


-- Only run if not called from another script
if (not pcall(debug.getlocal, 4, 1)) then
  unsatSolve(ARGS.problemID, ARGS.numPreRefs, ARGS.numRefs, ARGS.adaptive)
end