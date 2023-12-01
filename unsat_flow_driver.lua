-- Lua-Script to perform an unsaturated density flow problem
-- Author: Niklas Conen

local myPath = ug_get_current_path()
package.path = package.path..";".. myPath.."/config/?.lua;".. myPath.."/?.lua;"..myPath.."/config/trench/?.lua;"

ug_load_script("../scripts/ug_util.lua")
ug_load_script("/unsat_flow_util.lua")
ug_load_script("../scripts/util/solver_util.lua")

util.CheckAndPrintHelp("unsaturated density flow problem");

-- Parameters every problem uses
-- problem specific parameters are in the config file
ARGS =
{
    problemID         = util.GetParam("--problem-id", "trench2D"),
    numPreRefs        = util.GetParamNumber("--numPreRefs", 1, "number of refinements before parallel distribution"),
    numRefs           = util.GetParamNumber("--numRefs", 4, "number of refinements after parallel distribution"),
    dt			      = util.GetParamNumber("--dt", 0.001), -- time step length
    newton            = util.HasParamOption("--newton", false),
}

local problem = require(ARGS.problemID)

InitUG(problem.domain.dim, AlgebraType("CPU", 1))

local dom = util.CreateAndDistributeDomain(problem.domain.grid, ARGS.numRefs, ARGS.numPreRefs,  {})

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
local startTime = problem.time.start
local endTime = problem.time.stop
local dt   = problem.time.dt
local dtMin = problem.time.dtmin
local dtMax = problem.time.dtmax
local TOL = problem.time.tol
local dtred = problem.time.dtred

if ARGS.newton then
    util.SolveNonlinearTimeProblem(disc.u, domainDisc, solver, disc.vtk, ARGS.problemID.."_",
    "ImplEuler", 1.0, startTime, endTime, dt, dtMin, dtred)
else
    -- LIMEX time-stepping

    -- Solvers config.
    local limexLSolver = {}
    local limexNLSolver = {}

    local limexConvCheck=ConvCheck(1, 1e-12, 1e-10, true)
    limexConvCheck:set_supress_unsuccessful(true)

    --local lsolveCheck = ConvCheck(1000, 1e-12, 1e-10)
    local nstages = 2

    for i=1,nstages do
        limexLSolver[i] = util.solver.CreateSolver(problem.linSolver)
        print(limexLSolver[i]:config_string())
        --limexLSolver[i]:set_convergence_check(lsolveCheck)

        limexNLSolver[i] = NewtonSolver()
        limexNLSolver[i]:set_linear_solver(limexLSolver[i])
        limexNLSolver[i]:set_convergence_check(limexConvCheck)
    end

    -- Setup for time integrator
    local limex = LimexTimeIntegrator(nstages)
    for i=1,nstages do
        limex:add_stage(i, limexNLSolver[i], domainDisc )
    end

    local weightedMetricSpace=CompositeSpace()
    --local spaceP = VelEnergyComponentSpace("p", 2, inst.coef.EnergyTensorFlow)
    --local spaceC = L2ComponentSpace("c", 2, inst.coef.Conductivity2)
    local spaceC = L2ComponentSpace("c", 2)
    --local spaceP = VelEnergyComponentSpace("p", 2, ConstUserMatrix(1.0))
    local spaceP = H1ComponentSpace("p", 2)

    weightedMetricSpace:add(spaceP)
    weightedMetricSpace:add(spaceC)

    local concErrorEst = CompositeGridFunctionEstimator()
    -- concErrorEst:add(weightedMetricSpace)
    concErrorEst:add(spaceP)
    concErrorEst:add(spaceC)

    limex:add_error_estimator(concErrorEst)
    limex:set_tolerance(TOL)
    limex:set_stepsize_safety_factor(0.8)
    limex:set_time_step(dt)
    limex:set_dt_min(dtMin)
    limex:set_dt_max(dtMax)
    limex:set_increase_factor(2.0)
    limex:enable_matrix_cache()
    --limex:disable_matrix_cache()
    limex:set_conservative(true)

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
    local filename = (problem.output.filename) or (ARGS.problemID)
    local vtkobserver = VTKOutputObserver(problem.output.file..filename..".vtk", disc.vtk)
    limex:attach_observer(vtkobserver)

    -- post process for saving time step size
    local luaObserver = LuaCallbackObserver()
    function luaPostProcess(step, time, currdt)
        print("")
        print("Integral over salt mass fraction: "..Integral(disc.CompositeSaltMass, disc.u))
        print("")
        print("Integral over fluid phase volume: "..Integral(disc.CompositeSaturation, disc.u))
        print("")
        print(">>>> TimeStep: "..step..","..time..","..currdt.." <<<<")
        print("")
        return 0;
    end
    luaObserver:set_callback("luaPostProcess")
    limex:attach_observer(luaObserver)

    local sw = CuckooClock()
    sw:tic()
    -- Solve problem.
    limexConvCheck:set_minimum_defect(3e-10)
    limex:apply(disc.u, endTime, disc.u, 0.0)

    print ("CDELTA="..sw:toc())
end
