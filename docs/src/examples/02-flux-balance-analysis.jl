
# # Flux balance analysis (FBA)

# We will use [`flux_balance_analysis`](@ref) and several related functions to
# find the optimal flux distribution in the *E. coli* "core" model.

# If it is not already present, download the model and load the package:
import Downloads: download

!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

# next, load the necessary packages

import COBREXA as X
import AbstractFBCModels as A # for the accessors
import JSONFBCModels as J # for the model type
import Tulip as T # use any JuMP supported optimizer
import GLPK as G

model = A.load(J.JSONFBCModel, "e_coli_core.json") # load the model

# run FBA on the model using default settings

vt = X.flux_balance_analysis(model, T.Optimizer)

@test isapprox(vt.objective, 0.8739, atol = TEST_TOLERANCE) #src

# Alternatively, a constraint tree can be passed in as well

ctmodel = X.fbc_model_constraints(model)

# We can also pass some modifications to the optimizer
# Except for `X.silence`, all other optimizer modifications 
# are the same as those in JuMP.
vt = X.flux_balance_analysis(
    ctmodel,
    G.Optimizer;
    modifications = [
        X.silence
        X.set_objective_sense(X.J.MAX_SENSE) # JuMP is called J inside COBREXA
        X.set_optimizer(T.Optimizer) # change to Tulip from GLPK
        X.set_optimizer_attribute("IPM_IterationsLimit", 110) # Tulip specific setting
    ],
)

@test isapprox(vt.objective, 0.8739, atol = TEST_TOLERANCE) #src

# We can also modify the model. The most explicit way to do this is 
# to make a new constraint tree representation of the model.

import ConstraintTrees as C

fermentation = ctmodel.fluxes.EX_ac_e.value + ctmodel.fluxes.EX_etoh_e.value

forced_mixed_fermentation =
    ctmodel * :fermentation^C.Constraint(fermentation, (10.0, 1000.0)) # new modified model is created

vt = X.flux_balance_analysis(forced_mixed_fermentation, T.Optimizer; modifications = [X.silence])

@test isapprox(vt.objective, 0.6337, atol = TEST_TOLERANCE) #src

# Models that cannot be solved return `nothing`. In the example below, the 
# underlying model is modified.

ctmodel.fluxes.ATPM.bound = (1000.0, 10000.0) # TODO make mutable

vt = X.flux_balance_analysis(ctmodel, T.Optimizer; modifications = [X.silence])

@test isnothing(vt) #src

# Models can also be piped into the analysis functions

ctmodel.fluxes.ATPM.bound = (8.39, 10000.0) # revert
vt = ctmodel |> X.flux_balance_analysis(T.Optimizer; modifications = [X.silence])

@test isapprox(vt.objective, 0.8739, atol = TEST_TOLERANCE) #src
