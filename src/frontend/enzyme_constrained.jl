
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
$(TYPEDEF)

A simple struct storing information about the isozyme composition, including
subunit stoichiometry and turnover numbers.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Isozyme
    gene_product_stoichiometry::Dict{String,Float64}
    kcat_forward::Maybe{Float64} = nothing
    kcat_backward::Maybe{Float64} = nothing
end

export Isozyme

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
constructing it manually by using [`with_enzyme_constraints`](@ref).
"""
function enzyme_constrained_flux_balance_analysis(
    model::A.AbstractFBCModel;
    reaction_isozymes::Dict{String,Dict{String,Isozyme}},
    gene_product_molar_masses::Dict{String,Float64},
    capacity_limits::Vector{Tuple{String,Vector{String},Float64}},
    optimizer,
    settings = [],
)
    ct = fbc_model_constraints(model)

    # allocate variables for everything (nb. += doesn't associate right here)
    ct =
        ct +
        :fluxes_forward^unsigned_positive_contribution_variables(ct.fluxes) +
        :fluxes_reverse^unsigned_negative_contribution_variables(ct.fluxes) +
        :gene_product_amounts^gene_product_variables(model) +
        :isozyme_amounts^sum(
            Symbol(rid)^C.variables(keys = Symbol.(keys(is)), bounds = C.Between(0, Inf))
            for (rid, is) in reaction_isozymes
        )

    # connect all parts with constraints
    ct =
        ct *
        :directional_flux_balance^sign_split_constraints(
            ct.fluxes_forward,
            ct.fluxes_reverse,
            ct.fluxes,
        ) *
        :isozyme_flux_balance^isozyme_flux_constraints(
            ct.fluxes_forward,
            ct.fluxes_reverse,
            ct.isozyme_amounts,
            reaction_isozymes,
        ) *
        :gene_product_isozyme_balance^gene_product_isozyme_constraints(
            ct.isozyme_amounts,
            ct.gene_amounts,
            reaction_isozymes,
        ) *
        :gene_product_capacity_limits^C.ConstraintTree(
            Symbol(id) => C.Constraint(
                value = sum(
                    ct.gene_product_amounts[gp].value * gene_product_molar_masses[gp]
                    for gp in gps
                ),
                bound = C.Between(0, limit),
            ) for (id, gps, limit) in capacity_limits
        )

    optimized_constraints(ct; objective = ct.objective.value, optimizer, settings)
end

export enzyme_constrained_flux_balance_analysis
