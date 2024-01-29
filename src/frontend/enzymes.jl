
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

Uses [`enzyme_constraints`](@ref) internally.
"""
function enzyme_constrained_flux_balance_analysis(
    model::A.AbstractFBCModel;
    reaction_isozymes::Dict{String,Dict{String,Isozyme}},
    gene_product_molar_masses::Dict{String,Float64},
    capacity::Union{Vector{Tuple{String,Vector{String},Float64}},Float64},
    optimizer,
    settings = [],
)
    constraints = flux_balance_constraints(model)

    constraints += enzyme_variables(
        constraints;
        reaction_isozymes,
    )
    
    constraints *= enzyme_constraints(constraints;
        reaction_isozymes,
        gene_product_molar_masses,
        capacity,
    )
    
    optimized_constraints(
        constraints;
        objective = constraints.objective.value,
        optimizer,
        settings,
    )
end

export enzyme_constrained_flux_balance_analysis

"""
$(TYPEDSIGNATURES)

Run a basic simplified enzyme-constrained flux balance analysis on `model`.
Requires `reaction_isozymes`, which is a mapping of reaction identifiers to
[`Isozyme`](@ref) descriptions, and `gene_product_molar_masses` which is a
mapping of gene products to their molar masses. Internally, the cheapest and
fastest isozyme (minimum of MW/kcat for each isozyme) is used in the capacity
bound.

`capacity` may be a single number, which sets the limit of protein required for
all the fluxes in `reaction_isozymes`. Alternatively, `capacity` may be a vector
of identifier-fluxes-limit triples that make a constraint (identified by the
given identifier) that limits the enzyme requirements of the listed fluxes to
the given limit.

Uses [`simplified_enzyme_constraints`](@ref) internally.
"""
function simplified_enzyme_constrained_flux_balance_analysis(
    model::A.AbstractFBCModel;
    reaction_isozymes::Dict{String,Dict{String,Isozyme}},
    gene_product_molar_masses::Dict{String,Float64},
    capacity::Union{Vector{Tuple{String,Vector{String},Float64}},Float64},
    optimizer,
    settings = [],
)
    constraints = flux_balance_constraints(model)

    constraints += unidirectional_fluxes(constraints)
    
    constraints *= simplified_enzyme_constraints(
        constraints;
        reaction_isozymes,
        gene_product_molar_masses,
        capacity,
    )

    optimized_constraints(
        constraints;
        objective = constraints.objective.value,
        optimizer,
        settings,
    )
end

export simplified_enzyme_constrained_flux_balance_analysis
