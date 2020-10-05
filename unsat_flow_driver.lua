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

local import = require(ARGS.problemID)

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

-- CreateDomain will create a domain and solve it in a single process
-- CreateAndDistributeDomain will create a domain and distributes it to multiple
-- processes, numPreRefs will be performed before, numRefs after distribution
local dom = util.CreateAndDistributeDomain(problem.domain.grid, ARGS.numRefs, ARGS.numPreRefs,  {})

-- saves the refined grid
-- SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), "refined.ugx", 0.1)

-- create approximation space.
local approxSpace = util.unsat.CreateApproxSpace(problem, dom)

-- Index ordering
OrderCuthillMcKee(approxSpace, true);


-- Creating the Domain discretisation for the problem
util.unsat.CreateDomainDisc(problem, approxSpace)

-- Solver Config