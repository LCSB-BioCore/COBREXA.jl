
# Copyright (c) 2021-2024, University of Luxembourg
# Copyright (c) 2021-2024, Heinrich-Heine University Duesseldorf
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
$(TYPEDSIGNATURES)

Run a basic enzyme constrained flux balance analysis on `model`. The enzyme
model is parameterized by `reaction_isozymes`, which is a mapping of reaction
IDs (those used in the fluxes of the model) to named [`Isozyme`](@ref)s.
Additionally, `gene_molar_masses` and `capacity_limitations` should be supplied.
The latter is a vector of tuples, where each tuple represents a distinct bound
as `(bound_id, genes_in_bound, protein_mass_bound)`. Typically, `model` has
bounded exchange reactions, which are unnecessary in enzyme constrained models.
Unbound these reactions by listing their IDs in `unconstrain_reactions`, which
makes them reversible. Optimization `settings` are directly forwarded.

In the event that your model requires more complex build steps, consider
constructing it manually by using [`add_enzyme_constraints!`](@ref).
"""
function enzyme_constrained_flux_balance_analysis(
    model::A.AbstractFBCModel,
    reaction_isozymes::Dict{String,Dict{String,SimpleIsozyme}},
    gene_molar_masses::Dict{String,Float64},
    capacity_limitations::Vector{Tuple{String,Vector{String},Float64}};
    optimizer,
    unconstrain_reactions = String[],
    settings = [],
)
    m = fbc_model_constraints(model)

    # create enzyme variables
    m += :enzymes^enzyme_variables(model)

    # add enzyme equality constraints (stoichiometry)
    m = add_enzyme_constraints!(m, reaction_isozymes)

    # add capacity limitations
    for (id, gids, cap) in capacity_limitations
        m *= Symbol(id)^enzyme_capacity(m.enzymes, gene_molar_masses, gids, cap)
    end

    for rid in Symbol.(unconstrain_reactions)
        m.fluxes[rid].bound = C.Between(-1000.0, 1000.0)
    end

    optimized_constraints(m; objective = m.objective.value, optimizer, settings)
end

export enzyme_constrained_flux_balance_analysis
