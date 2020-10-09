-- Lua-Script to perform an unsaturated density flow problem
-- Author: Niklas Conen

local myPath = ug_get_current_path()
package.path = package.path..";".. myPath.."/config/?.lua;".. myPath.."/?.lua"

ug_load_script("../scripts/ug_util.lua")
ug_load_script("/unsat_flow_util.lua")


ARGS = {
  problemID = util.GetParam("--problem-id", "trench2D"),
  numPreRefs = util.GetParamNumber("--numPreRefs", 1, "number of refinements before parallel distribution"),
  numRefs = util.GetParamNumber("--num-refs", 2, "number of refinements after parallel distribution"),
  check = util.HasParamOption("--check", false, "checks if the config file has the correct layout")
  }

util.CheckAndPrintHelp("unsaturated density flow problem");

local problem = require(ARGS.problemID)

InitUG(problem.domain.dim, AlgebraType("CPU", 1))

local vtk = VTKOutput()

if ARGS.check then
    local check = require("problem_check")
    print("checking problem file ", ARGS.problemID, ".lua")
    if not problem_check(problem) then
        print("problem check failed.")
        exit()
    end
    print()
end

local dom = util.CreateAndDistributeDomain(problem.domain.grid, ARGS.numRefs, ARGS.numPreRefs,  {})

-- saves the refined grid
-- SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), "refined.ugx", 0.1)

-- create approximation space.
local approxSpace = util.unsat.CreateApproxSpace(problem, dom)

-- Index ordering
OrderCuthillMcKee(approxSpace, true);

-- Creating the Domain discretisation for the problem
local domainDisc = util.unsat.CreateDomainDisc(problem, approxSpace)

-- Solver Config

local u = GridFunction(approxSpace)
util.unsat.SetInitialData(problem, u)


local linSolver = LU()

local newtonConvCheck = ConvCheck()
newtonConvCheck:set_maximum_steps(10)
newtonConvCheck:set_minimum_defect(5e-8)
newtonConvCheck:set_reduction(1e-10)
newtonConvCheck:set_verbose(true)

local newtonLineSearch = StandardLineSearch()

local newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(linSolver)
newtonSolver:set_convergence_check(newtonConvCheck)
newtonSolver:set_line_search(newtonLineSearch)

-- Time stepping parameters
local endTime = problem.time.stop
local dt   = problem.time.dt
local dtMin = problem.time.dtmin
local dtMax = 1e-1*endTime  -- at least 10%
local TOL = problem.time.tol or ARGS.limexTOL

-- since solver configurations can be quite complex, we print the configuration:
print("NewtonSolver configuration:")
print(newtonSolver:config_string())

newtonConvCheck:set_maximum_steps(100)
util.SolveNonlinearTimeProblem(u, domainDisc, newtonSolver, vtk, "Newton",
"ImplEuler", 1.0, 0.0,  dt, dt, dtMin, 0.5)