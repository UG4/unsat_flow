"""
Applies the ParameterEstimator Plugin to estimate the hydraulic conductivity K, 
the porosity n_p, and van Genuchten Parameters alpha and n

Parameter Variations:
    - K: 2e-3 to 1e-2 m/s
    - n_p: 0.25 to 0.4
    - alpha: 1.3 to 20 1/m
    - n: 0.3 to 3
"""

# import from plugins folder
# please make sure the enivorment variable UG4_ROOT points to your UG4 directory!
import sys, os
sys.path.append(os.path.join(os.environ["UG4_ROOT"],"plugins","ParameterEstimator", "python")) 

from UGParameterEstimator import *

# specify the parameters
pm = ParameterManager()

# direct parameter: no transformation from parameter used in optimizer to parameter in lua
# the specified value is the initial value used in the optimization
# lower and upper bounds can be supplied as third and fourth argument
pm.addParameter(DirectParameter("thetaS", 0.35, 0.25, 0.4))
pm.addParameter(DirectParameter("Ksat", 2.2e-3, 2e-3, 1e-2))

# create the evaluator object. this will create a LocalEvaluator, which uses MPI for parallelism > 1, if no UGSUBMIT was found,
# or an ClusterEvaluator using UGSUBMIT.
evaluator = Evaluator.ConstructEvaluator(
    luafile="evaluate.lua",             # the lua file to execute for every evaluation
    directory="evaluations",            # the folder used for data exchange
    parametermanager=pm,                # the parameters defined above
    evaluation_type=GenericEvaluation,         # the type the evaluations should be parsed as.
    parameter_output_adapter=UG4ParameterOutputAdapter(),       # the adapter to use to write the parameters
    threadcount=1)                     # threads to use locally or in ugsubmit when using UGSUBMIT

# create the optimizer
optimizer = GaussNewtonOptimizer(LogarithmicParallelLineSearch(evaluator))

# specify some fixed parameters if needed (could be done in lua, also)
evaluator.fixedparameters["n"] = 2.06
evaluator.fixedparameters["alpha"] = 0.423/(998.23*9.81)

# this will do a measurement with fixed parameters
target = GenericEvaluation([], []).fromCSV("target.csv")

# try to restore these parameters by calibration
# store the calibration process and logging in example.pkl
result = Result("pump.pkl")
result = optimizer.run(evaluator, pm.getInitialArray(), target, result=result)
