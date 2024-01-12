
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

# # Loopless flux balance analysis (ll-FBA)

# Here we wil add loopless constraints to a flux balance model to ensure that
# the resultant solution is thermodynamically consistent. As before, we will use
# the core *E. coli* model, which we can download using
# [`download_model`](@ref):

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

# Additionally to COBREXA and the JSON model format package. We will also need a
# solver that can solve mixed interger linear programs like GLPK.

import JSONFBCModels
import GLPK
import AbstractFBCModels as A

model = load_model("e_coli_core.json")

# ## Running a simple loopless FBA (ll-FBA)

# One can directly use `loopless_flux_balance_analysis` to solve an FBA problem
# based on `model` where loopless constraints are added to all fluxes. This is
# the direct approach.

sol = loopless_flux_balance_analysis(model; optimizer = GLPK.Optimizer)

@test isapprox(sol.objective, 0.8739215069684303, atol = TEST_TOLERANCE) #src

@test all(
    v * sol.pseudo_gibbs_free_energy_reaction[k] <= -TEST_TOLERANCE for
    (k, v) in sol.fluxes if
    haskey(sol.pseudo_gibbs_free_energy_reaction, k) && abs(v) >= 1e-6
) #src

# ## Building your own loopless model

# ConstraintTrees allows one to add loopless constraints to any model. To
# illustrate how one would add loopless constraints to an arbitrary model (and
# not use the convenience function), let's build a loopless model from scratch.

# First, build a normal flux balance model
m = flux_balance_constraints(model)

# Next, find all internal reactions, and their associated indices for use later
internal_reactions = [
    (i, Symbol(rid)) for
    (i, rid) in enumerate(A.reactions(model)) if !is_boundary(model, rid)
]
internal_reaction_ids = last.(internal_reactions)
internal_reaction_idxs = first.(internal_reactions) # order needs to match the internal reaction ids below

# Construct the stoichiometric nullspace of the internal reactions
import LinearAlgebra: nullspace

internal_reaction_stoichiometry_nullspace_columns =
    eachcol(nullspace(Array(A.stoichiometry(model)[:, internal_reaction_idxs])))

# And simply add loopless contraints on the fluxes of the model
m = with_loopless_constraints(
    m,
    internal_reaction_ids,
    internal_reaction_stoichiometry_nullspace_columns;
    fluxes = m.fluxes,
)

# Now the model can be solved as before!
optimized_constraints(m; objective = m.objective.value, optimizer = GLPK.Optimizer)
