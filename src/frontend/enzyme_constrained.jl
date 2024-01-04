
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
    kcat_reverse::Maybe{Float64} = nothing
end

export Isozyme

"""
$(TYPEDSIGNATURES)

Run a basic enzyme-constrained flux balance analysis on `model`. The enzyme
model is parameterized by `reaction_isozymes`, which is a mapping of reaction
identifiers to [`Isozyme`](@ref) descriptions.

Additionally, one typically wants to supply `gene_product_molar_masses` to
describe the weights of enzyme building material, and `capacity` which limits
the mass of enzymes in the whole model.

`capacity` may be a single number, which sets the limit for "all described
enzymes". Alternatively, `capacity` may be a vector of identifier-genes-limit
triples that make a constraint (identified by the given identifier) that limits
the listed genes to the given limit.
"""
function enzyme_constrained_flux_balance_analysis(
    model::A.AbstractFBCModel;
    reaction_isozymes::Dict{String,Dict{String,Isozyme}},
    gene_product_molar_masses::Dict{String,Float64},
    capacity::Union{Vector{Tuple{String,Vector{String},Float64}},Float64},
    optimizer,
    settings = [],
)
    ct = fbc_model_constraints(model)

    # might be nice to omit some conditionally (e.g. slash the direction if one
    # kcat is nothing)
    isozyme_amounts = isozyme_amount_variables(
        Symbol.(keys(reaction_isozymes)),
        rid -> Symbol.(keys(reaction_isozymes[string(rid)])),
    )

    # allocate variables for everything (nb. += wouldn't associate right here)
    ct =
        ct +
        :fluxes_forward^unsigned_positive_contribution_variables(ct.fluxes) +
        :fluxes_reverse^unsigned_negative_contribution_variables(ct.fluxes) +
        :isozyme_forward_amounts^isozyme_amounts +
        :isozyme_reverse_amounts^isozyme_amounts +
        :gene_product_amounts^C.variables(
            keys = Symbol.(A.genes(model)),
            bounds = Between(0, Inf),
        )

    # connect all parts with constraints
    ct =
        ct *
        :directional_flux_balance^sign_split_constraints(
            ct.fluxes_forward,
            ct.fluxes_reverse,
            ct.fluxes,
        ) *
        :isozyme_flux_forward_balance^isozyme_flux_constraints(
            ct.isozyme_forward_amounts,
            ct.fluxes_forward,
            (rid, isozyme) -> maybemap(
                x -> x.kcat_forward,
                maybeget(reaction_isozymes, string(rid), string(isozyme)),
            ),
        ) *
        :isozyme_flux_reverse_balance^isozyme_flux_constraints(
            ct.isozyme_reverse_amounts,
            ct.fluxes_reverse,
            (rid, isozyme) -> maybemap(
                x -> x.kcat_reverse,
                maybeget(reaction_isozymes, string(rid), string(isozyme)),
            ),
        ) *
        :gene_product_isozyme_balance^gene_product_isozyme_constraints(
            ct.gene_product_amounts,
            (ct.isozyme_forward_amounts, ct.isozyme_reverse_amounts),
            (rid, isozyme) -> maybemap(
                x -> x.gene_product_stoichiometry,
                maybeget(reaction_isozymes, string(rid), string(isozyme)),
            ),
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
