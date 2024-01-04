
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

Allocate a (non-negative) variable for all amounts of gene products in the
`model`.
"""
gene_product_variables(model::A.AbstractFBCModel) =
    C.variables(; keys = Symbol.(A.genes(model)), bounds = C.Between(0.0, Inf))

export gene_product_variables

"""
$(TYPEDSIGNATURES)

Helper function to link isozymes to fluxes. All forward (backward) fluxes are
the sum of all the isozymes catalysing these fluxes.
"""
function link_isozymes(
    fluxes_directional::C.ConstraintTree,
    fluxes_isozymes::C.ConstraintTree,
)
    C.ConstraintTree(
        k => C.Constraint(
            value = s.value - sum(x.value for (_, x) in fluxes_isozymes[k]),
            bound = C.EqualTo(0.0),
        ) for (k, s) in fluxes_directional if haskey(fluxes_isozymes, k)
    )
end

"""
$(TYPEDSIGNATURES)

Helper function to create the enzyme "mass balance" matrix. In essence, the
stoichiometric coefficient is subunit_stoichiometry / kcat for each directional,
isozyme flux, and it must be balanced by the enzyme variable supply.
"""
function enzyme_stoichiometry(
    enzymes::C.ConstraintTree,
    fluxes_isozymes_forward::C.ConstraintTree,
    fluxes_isozymes_backward::C.ConstraintTree,
    reaction_isozymes::Dict{String,Dict{String,T}},
) where {T<:Isozyme}
    # map enzyme ids to reactions that use them (not all isozymes have to though)
    enzyme_rid_lookup = Dict{Symbol,Vector{Symbol}}()
    for (rid, isozymes) in reaction_isozymes
        for isozyme in values(isozymes)
            for gid in keys(isozyme.gene_product_stoichiometry)
                rids = get!(enzyme_rid_lookup, Symbol(gid), Symbol[])
                rid in rids || push!(rids, Symbol(rid))
            end
        end
    end

    C.ConstraintTree(
        gid => C.Constraint(
            value = enz.value + # supply
                    sum(
                        enzyme_balance(
                            gid,
                            rid,
                            fluxes_isozymes_forward,
                            reaction_isozymes,
                            :kcat_forward,
                        ) for rid in enzyme_rid_lookup[gid] if
                        haskey(fluxes_isozymes_forward, rid);
                        init = zero(typeof(enz.value)),
                    ) + # flux through positive isozymes
                    sum(
                        enzyme_balance(
                            gid,
                            rid,
                            fluxes_isozymes_backward,
                            reaction_isozymes,
                            :kcat_backward,
                        ) for rid in enzyme_rid_lookup[gid] if
                        haskey(fluxes_isozymes_backward, rid);
                        init = zero(typeof(enz.value)),
                    ), # flux through negative isozymes
            bound = C.EqualTo(0.0),
        ) for (gid, enz) in enzymes if gid in keys(enzyme_rid_lookup)
    )
end

"""
$(TYPEDSIGNATURES)

Helper function to balance the forward or backward isozyme fluxes for a specific
gene product.
"""
function enzyme_balance(
    gid::Symbol,
    rid::Symbol,
    fluxes_isozymes::C.ConstraintTree, # direction
    reaction_isozymes::Dict{String,Dict{String,T}},
    direction = :kcat_forward,
) where {T<:Isozyme}
    isozyme_dict = Dict(Symbol(k) => v for (k, v) in reaction_isozymes[string(rid)])

    sum( # this is where the stoichiometry comes in
        -isozyme_value.value *
        isozyme_dict[isozyme_id].gene_product_stoichiometry[string(gid)] /
        getfield(isozyme_dict[isozyme_id], direction) for
        (isozyme_id, isozyme_value) in fluxes_isozymes[rid] if
        gid in Symbol.(keys(isozyme_dict[isozyme_id].gene_product_stoichiometry));
        init = zero(C.LinearValue),
    )
end

"""
$(TYPEDSIGNATURES)

Create a enzyme capacity limitation. Bounds the gene product masses (concentration
* molar mass) of gene products in `enzyme_ids` by `capacity_bound`.
"""
function enzyme_capacity(
    enzymes::C.ConstraintTree,
    gene_molar_masses::Dict{String,Float64},
    enzyme_ids::Vector{String},
    capacity_bound::C.Bound,
)
    C.Constraint(
        value = sum(
            enzymes[Symbol(gid)].value * gene_molar_masses[gid] for gid in enzyme_ids
        ),
        bound = capacity_bound,
    )
end

"""
$(TYPEDSIGNATURES)

Create an enzyme capacity limitation. Bounds the gene product masses
(concentration * molar mass) of gene products in `enzyme_ids` between `[0, capacity]`.
"""
enzyme_capacity(
    enzymes::C.ConstraintTree,
    gene_molar_masses::Dict{String,Float64},
    enzyme_ids::Vector{String},
    capacity::Float64,
) = enzyme_capacity(enzymes, gene_molar_masses, enzyme_ids, C.Between(0.0, capacity))


export enzyme_capacity

"""
$(TYPEDSIGNATURES)

Add enzyme constraints to a constraint tree, `m`. The enzyme model is
parameterized by `reaction_isozymes`, which is a mapping of reaction IDs (those
used in the fluxes of the model) to named struct which is a subtype of
[`Isozyme`](@ref)s. Additionally, `gene_molar_masses` and `capacity_limitations`
should be supplied. The latter is a vector of tuples, where each tuple
represents a distinct bound as `(bound_id, genes_in_bound, protein_mass_bound)`.
Finally, specify the `fluxes` and `enzymes` to which the constraints should be
mounted.

# Note
The isozyme struct used in `reaction_isozymes` must have fields
`gene_product_stoichiometry`, `kcat_forward`, and `kcat_backward` to properly
assign kcats to reactions. Use [`Isozyme`](@ref) when in doubt.
"""
function with_enzyme_constraints(
    m::C.ConstraintTree,
    reaction_isozymes::Dict{String,Dict{String,Isozyme}};
    fluxes = m.fluxes,
    enzymes = m.enzymes,
)

    # create directional fluxes
    m +=
        :fluxes_forward^fluxes_in_direction(fluxes, :forward) +
        :fluxes_backward^fluxes_in_direction(fluxes, :backward)

    # link directional fluxes to original fluxes
    m *=
        :link_flux_directions^sign_split_constraints(
            positive = m.fluxes_forward,
            negative = m.fluxes_backward,
            signed = fluxes,
        )

    # create fluxes for each isozyme
    for (rid, _) in m.fluxes_forward
        if haskey(reaction_isozymes, string(rid))
            m +=
                :fluxes_isozymes_forward^rid^C.variables(
                    keys = Symbol.(keys(reaction_isozymes[string(rid)])),
                    bounds = Between(0, Inf),
                )
        end
    end
    for (rid, _) in m.fluxes_backward
        if haskey(reaction_isozymes, string(rid))
            m +=
                :fluxes_isozymes_backward^rid^C.variables(
                    keys = Symbol.(keys(reaction_isozymes[string(rid)])),
                    bounds = Between(0, Inf),
                )
        end
    end

    # link isozyme fluxes to directional fluxes
    m *=
        :link_isozyme_fluxes_forward^link_isozymes(
            m.fluxes_forward,
            m.fluxes_isozymes_forward,
        )
    m *=
        :link_isozyme_fluxes_backward^link_isozymes(
            m.fluxes_backward,
            m.fluxes_isozymes_backward,
        )

    # add enzyme mass balances
    m *=
        :enzyme_stoichiometry^enzyme_stoichiometry(
            enzymes,
            m.fluxes_isozymes_forward,
            m.fluxes_isozymes_backward,
            reaction_isozymes,
        )

    m
end

export with_enzyme_constraints
