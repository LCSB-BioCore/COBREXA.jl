
# # Parsimonious flux balance analysis

# We will use [`parsimonious_flux_balance`](@ref) and
# [`minimize_metabolic_adjustment`](@ref) to find the optimal flux
# distribution in the *E. coli* "core" model.
#
# TODO pFBA citation

# If it is not already present, download the model and load the package:
import Downloads: download

!isfile("e_coli_core.json") &&
    download("http://bigg.ucsd.edu/static/models/e_coli_core.json", "e_coli_core.json")

# next, load the necessary packages

using COBREXA

import JSONFBCModels
import Clarabel # can solve QPs

model = load_model("e_coli_core.json") # load the model

# Use the convenience function to run standard pFBA on

vt =
    parsimonious_flux_balance_analysis(model, Clarabel.Optimizer; modifications = [silence])

# Or use the piping functionality

model |> parsimonious_flux_balance_analysis(Clarabel.Optimizer; modifications = [silence])

@test isapprox(vt.objective, 0.87392; atol = TEST_TOLERANCE) #src
@test sum(x^2 for x in values(vt.fluxes)) < 15000 #src

#=

# Alternatively, you can construct your own constraint tree model with
# the quadratic objective (this approach is much more flexible).

ctmodel = build_flux_balance_model(model)
ctmodel *= :l2objective^squared_sum_objective(ctmodel.fluxes)
ctmodel.objective.bound = 0.3 # set growth rate # TODO currently breaks

opt_model = optimization_model(
    ctmodel;
    objective = ctmodel.:l2objective.value,
    optimizer = Clarabel.Optimizer,
    sense = Minimal,
)

J.optimize!(opt_model) # JuMP is called J in COBREXA

is_solved(opt_model) # check if solved

vt = C.constraint_values(ctmodel, J.value.(opt_model[:x])) # ConstraintTrees.jl is called C in COBREXA

@test isapprox(vt.l2objective, ?; atol = QP_TEST_TOLERANCE) #src  # TODO will break until mutable bounds

# It is likewise as simple to run MOMA using the convenience functions.

ref_sol = Dict("ATPS4r" => 33.0, "CYTBD" => 22.0)

vt = minimize_metabolic_adjustment(model, ref_sol, Gurobi.Optimizer)

# Or use the piping functionality

model |>
minimize_metabolic_adjustment(ref_sol, Clarabel.Optimizer; modifications = [silence])

@test isapprox(vt.:momaobjective, 0.81580806; atol = TEST_TOLERANCE) #src

# Alternatively, you can construct your own constraint tree model with
# the quadratic objective (this approach is much more flexible).

ctmodel = build_flux_balance_model(model)
ctmodel *=
    :minoxphospho^squared_sum_error_objective(
        ctmodel.fluxes,
        Dict(:ATPS4r => 33.0, :CYTBD => 22.0),
    )
ctmodel.objective.bound = 0.3 # set growth rate # TODO currently breaks

opt_model = optimization_model(
    ctmodel;
    objective = ctmodel.minoxphospho.value,
    optimizer = Clarabel.Optimizer,
    sense = Minimal,
)

J.optimize!(opt_model) # JuMP is called J in COBREXA

is_solved(opt_model) # check if solved

vt = C.constraint_values(ctmodel, J.value.(opt_model[:x])) # ConstraintTrees.jl is called C in COBREXA

@test isapprox(vt.l2objective, ?; atol = QP_TEST_TOLERANCE) #src  # TODO will break until mutable bounds

=#
