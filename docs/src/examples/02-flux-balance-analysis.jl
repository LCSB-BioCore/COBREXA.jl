
# # Flux balance analysis

import COBREXA as X
import AbstractFBCModels as A
import JSONFBCModels as J

model = A.load(J.JSONFBCModel, "e_coli_core.json")

ctmodel = X.fbc_model_constraints(model)


