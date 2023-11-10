
# # Flux balance analysis

import COBREXA as X
import AbstractFBCModels as A
import JSONFBCModels as J
import ConstraintTrees as C
using Gurobi

model = A.load(J.JSONFBCModel, "e_coli_core.json")

ctmodel = X.fbc_model_constraints(model)

vt = C.ValueTree(ctmodel, value.(opt_model[:x]))

vt = X.flux_balance_analysis(model, Gurobi.Optimizer)

vt = X.flux_balance_analysis(ctmodel, Gurobi.Optimizer; modifications=[X.silence])
