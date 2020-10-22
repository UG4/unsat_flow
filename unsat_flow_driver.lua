-- Lua-Script to perform an unsaturated density flow problem
-- Author: Niklas Conen

local myPath = ug_get_current_path()
package.path = package.path..";".. myPath.."/config/?.lua;".. myPath.."/?.lua"

ug_load_script("../scripts/ug_util.lua")
ug_load_script("/unsat_flow_util.lua")


ARGS = {
  problemID = util.GetParam("--problem-id", "trench2D"),
  numPreRefs = util.GetParamNumber("--numPreRefs", 1, "number of refinements before parallel distribution"),
  numRefs = util.GetParamNumber("--numRefs", 2, "number of refinements after parallel distribution"),
  check = util.HasParamOption("--check", false, "checks if the config file has the correct layout"),
  outFileNamePrefix = util.GetParam("-o", "unsat_"),
  }

util.CheckAndPrintHelp("unsaturated density flow problem");

local problem = require(ARGS.problemID)

if ARGS.check then
    local check = require("problem_check")
    print("checking problem file ", ARGS.problemID, ".lua")
    if not problem_check(problem) then
        print("problem check failed.")
        exit()
    end
    print()
end


InitUG(problem.domain.dim, AlgebraType("CPU", 1))

local vtk = VTKOutput()

local dom = util.CreateAndDistributeDomain(problem.domain.grid, ARGS.numRefs, ARGS.numPreRefs,  {})

-- saves the refined grid
-- SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), "refined.ugx", 0.1)

disc = ProblemDisc:new(problem, dom, vtk)

-- create approximation space.
approxSpace = disc:CreateApproxSpace()

-- Index ordering
OrderCuthillMcKee(approxSpace, true);

-- Creating the Domain discretisation for the problem
domainDisc = disc:CreateDomainDisc(approxSpace)

-- Initial Data
u = GridFunction(approxSpace)
disc:SetInitialData(u)

-- Solver Config
solver = disc:CreateSolver()

-- Time stepping parameters
local endTime = problem.time.stop
local dt   = problem.time.dt
local dtMin = problem.time.dtmin
local dtMax = 1e-1*endTime  -- at least 10%
local TOL = problem.time.tol or ARGS.limexTOL

util.SolveNonlinearTimeProblem(u, domainDisc, solver, vtk, ARGS.outFileNamePrefix,
"ImplEuler", 1.0, 0.0,  endTime, dt, dtMin, 0.5)

vtk:print(filename, u, step, time)
