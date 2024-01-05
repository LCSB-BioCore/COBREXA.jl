
# Copyright (c) 2021-2024, University of Luxembourg                         #src
# Copyright (c) 2021-2024, Heinrich-Heine University Duesseldorf            #src
#                                                                           #src
# Licensed under the Apache License, Version 2.0 (the "License");           #src
# you may not use this file except in compliance with the License.          #src
# You may obtain a copy of the License at                                   #src
#                                                                           #src
#     http://www.apache.org/licenses/LICENSE-2.0                            #src
#                                                                           #src
# Unless required by applicable law or agreed to in writing, software       #src
# distributed under the License is distributed on an "AS IS" BASIS,         #src
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  #src
# See the License for the specific language governing permissions and       #src
# limitations under the License.                                            #src

# # Parsimonious flux balance analysis

# We will use [`parsimonious_flux_balance_analysis`](@ref) and
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

vt = parsimonious_flux_balance_analysis(model, Clarabel.Optimizer; settings = [silence])

# Or use the piping functionality

model |> parsimonious_flux_balance_analysis(Clarabel.Optimizer; settings = [silence])

@test isapprox(vt.objective, 0.87392; atol = TEST_TOLERANCE) #src
@test sum(x^2 for x in values(vt.fluxes)) < 15000 #src

#=

# Alternatively, you can construct your own constraint tree model with
# the quadratic objective (this approach is much more flexible).

ctmodel = fbc_model_constraints(model)
ctmodel *= :l2objective^squared_sum_value(ctmodel.fluxes)
ctmodel.objective.bound = 0.3 # set growth rate # TODO currently breaks

opt_model = optimization_model(
    ctmodel;
    objective = ctmodel.:l2objective.value,
    optimizer = Clarabel.Optimizer,
    sense = Minimal,
)

J.optimize!(opt_model) # JuMP is called J in COBREXA

is_solved(opt_model) # check if solved

vt = C.substitute_values(ctmodel, J.value.(opt_model[:x])) # ConstraintTrees.jl is called C in COBREXA

@test isapprox(vt.l2objective, ?; atol = QP_TEST_TOLERANCE) #src  # TODO will break until mutable bounds

# It is likewise as simple to run MOMA using the convenience functions.

ref_sol = Dict("ATPS4r" => 33.0, "CYTBD" => 22.0)

vt = minimize_metabolic_adjustment(model, ref_sol, Gurobi.Optimizer)

# Or use the piping functionality

model |>
minimize_metabolic_adjustment(ref_sol, Clarabel.Optimizer; settings = [silence])

@test isapprox(vt.:momaobjective, 0.81580806; atol = TEST_TOLERANCE) #src

# Alternatively, you can construct your own constraint tree model with
# the quadratic objective (this approach is much more flexible).

ctmodel = fbc_model_constraints(model)
ctmodel *=
    :minoxphospho^squared_sum_error_value(
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

vt = C.substitute_values(ctmodel, J.value.(opt_model[:x])) # ConstraintTrees.jl is called C in COBREXA

@test isapprox(vt.l2objective, ?; atol = QP_TEST_TOLERANCE) #src  # TODO will break until mutable bounds

=#
