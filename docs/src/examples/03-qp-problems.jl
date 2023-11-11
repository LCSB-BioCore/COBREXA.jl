
# # Quandratic objective flux balance analysis type problems

# We will use [`parsimonious_flux_balance_analysis`](@ref) and
# [`minimize_metabolic_adjustment_analysis`](@ref) to find the optimal flux
# distribution in the *E. coli* "core" model.

# If it is not already present, download the model and load the package:
import Downloads: download

!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

# next, load the necessary packages

import COBREXA as X
import AbstractFBCModels as A # for the accessors
import JSONFBCModels as J # for the model type
import Clarabel # can solve QPs

model = A.load(J.JSONFBCModel, "e_coli_core.json") # load the model

# Use the convenience function to run standard pFBA

vt = X.parsimonious_flux_balance_analysis(model, Clarabel.Optimizer)

# Or use the piping functionality

model |> parsimonious_flux_balance_analysis(Clarabel.Optimizer; modifications = [X.silence])

@test isapprox(vt.objective, 0.87392; atol=TEST_TOLERANCE) #src
@test isapprox(sum(x^2 for x in values(vt.fluxes)), 11414.2119; atol=QP_TEST_TOLERANCE) #src

# Alternatively, you can construct your own constraint tree model with 
# the quadratic objective (this approach is much more flexible).

ctmodel = X.fbc_model_constraints(model)
ctmodel *= :l2objective ^ X.squared_sum_objective(ctmodel.fluxes)
ctmodel.objective.bound = 0.3 # set growth rate # TODO currently breaks

opt_model = X.optimization_model(
    ctmodel;
    objective = ctmodel.:l2objective.value,
    optimizer = Clarabel.Optimizer,
    sense = X.J.MIN_SENSE,
)

X.J.optimize!(opt_model) # JuMP is called J in COBREXA
    
X.is_solved(opt_model) # check if solved

vt = X.C.ValueTree(ctmodel, X.J.value.(opt_model[:x])) # ConstraintTrees.jl is called C in COBREXA

@test isapprox(vt.l2objective, ?; atol=QP_TEST_TOLERANCE) #src  # TODO will break until mutable bounds

# It is likewise as simple to run MOMA using the convenience functions. 

ref_sol = Dict("ATPS4r" => 33.0, "CYTBD" => 22.0)

vt = X.minimize_metabolic_adjustment_analysis(model, ref_sol, Gurobi.Optimizer)

# Or use the piping functionality

model |> X.minimize_metabolic_adjustment_analysis(ref_sol, Clarabel.Optimizer; modifications = [X.silence])

@test isapprox(vt.:momaobjective, 0.81580806; atol=TEST_TOLERANCE) #src

# Alternatively, you can construct your own constraint tree model with 
# the quadratic objective (this approach is much more flexible).

ctmodel = X.fbc_model_constraints(model)
ctmodel *= :minoxphospho ^ X.squared_sum_error_objective(ctmodel.fluxes, Dict(:ATPS4r => 33.0, :CYTBD => 22.0,))
ctmodel.objective.bound = 0.3 # set growth rate # TODO currently breaks

opt_model = X.optimization_model(
    ctmodel;
    objective = ctmodel.minoxphospho.value,
    optimizer = Clarabel.Optimizer,
    sense = X.J.MIN_SENSE,
)

X.J.optimize!(opt_model) # JuMP is called J in COBREXA
    
X.is_solved(opt_model) # check if solved

vt = X.C.ValueTree(ctmodel, X.J.value.(opt_model[:x])) # ConstraintTrees.jl is called C in COBREXA

@test isapprox(vt.l2objective, ?; atol=QP_TEST_TOLERANCE) #src  # TODO will break until mutable bounds
