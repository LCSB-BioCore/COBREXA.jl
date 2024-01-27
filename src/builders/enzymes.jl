
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
subunit stoichiometry and turnover numbers. Use with
[`enzyme_constrained_flux_balance_analysis`](@ref).

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

Create a `ConstraintTree` with variables for isozyme contributions to reaction
fluxes. The tree has 2 levels: the first contains all reaction flux IDs that
have isozymes, the second contains the isozyme IDs for each reaction flux.

`fluxes` should be anything that can be iterated to give reaction flux IDs.

`flux_isozymes` is a function that, for a given reaction flux ID, returns
anything iterable that contains the isozyme IDs for the given reaction flux.
Returning an empty iterable prevents allocating the subtree for the given flux.
"""
isozyme_amount_variables(fluxes, flux_isozymes) = sum(
    (
        f^C.variables(keys = flux_isozymes(f), bounds = C.Between(0, Inf)) for
        f in fluxes if !isempty(flux_isozymes(f))
    ),
    init = C.ConstraintTree(),
)

export isozyme_amount_variables

"""
$(TYPEDSIGNATURES)

A constraint tree that sums up partial contributions of reaction isozymes to
the fluxes of reactions.

For practical purposes, both fluxes and isozymes are here considered to be
unidirectional, i.e., one would typically apply this twice to constraint both
"forward" and "reverse" fluxes.

Function `kcat` should retugn the kcat value for a given reaction and isozyme
(IDs of which respectively form the 2 parameters for each call).
"""
function isozyme_flux_constraints(
    isozyme_amounts::C.ConstraintTree,
    fluxes::C.ConstraintTree,
    kcat,
)
    C.ConstraintTree(
        rid => C.Constraint(
            sum(kcat(rid, iid) * i.value for (iid, i) in ri if !isnothing(kcat(rid, iid))) - fluxes[rid].value,
            0.0,
        ) for (rid, ri) in isozyme_amounts if haskey(fluxes, rid)
    )
end

export isozyme_flux_constraints

"""
$(TYPEDSIGNATURES)

A constraint tree that binds the isozyme amounts to gene product amounts
accordingly to their multiplicities (aka. stoichiometries, protein units, ...)
given by `isozyme_stoichiometry`.

Values in `gene_product_amounts` should describe the gene product allocations.

`isozyme_amount_trees` is an iterable that contains `ConstraintTree`s that
describe the allocated isozyme amounts (such as created by
[`isozyme_amount_variables`](@ref). One may typically pass in both forward- and
reverse-direction amounts at once, but it is possible to use a single tree,
e.g., in a uni-tuple: `tuple(my_tree)`.

`isozyme_stoichiometry` gets called with a reaction and isozyme ID as given by
the isozyme amount trees. It may return `nothing` in case there's no
information.
"""
function gene_product_isozyme_constraints(
    gene_product_amounts::C.ConstraintTree,
    isozymes_amount_trees,
    isozyme_stoichiometry,
)
    res = C.ConstraintTree()
    # This needs to invert the stoichiometry mapping,
    # so we patch up a fresh new CT in place.
    for iss in isozymes_amount_trees
        for (rid, is) in iss
            for (iid, i) in is
                gpstoi = isozyme_stoichiometry(rid, iid)
                isnothing(gpstoi) && continue
                for (gp, stoi) in gpstoi
                    haskey(gene_product_amounts, gp) || continue
                    if haskey(res, gp)
                        res[gp].value += i.value * stoi
                    else
                        res[gp] =
                            C.Constraint(i.value * stoi - gene_product_amounts[gp].value, 0)
                    end
                end
            end
        end
    end
    res
end

export gene_product_isozyme_constraints

function add_enzyme_constraints(
    constraints::C.ConstraintTree;
    fluxes = constraints.fluxes,
    reaction_isozymes::Dict{String,Dict{String,Isozyme}},
    gene_product_molar_masses::Dict{String,Float64},
    capacity::Union{Vector{Tuple{String,Vector{String},Float64}},Float64},
)

    gene_ids = unique([
        Symbol(gid) for isos_dict in values(reaction_isozymes) for
        iso in values(isos_dict) for gid in keys(iso.gene_product_stoichiometry)
    ])

    # might be nice to omit some conditionally (e.g. slash the direction if one
    # kcat is nothing)
    isozyme_amounts = isozyme_amount_variables(
        Symbol.(keys(reaction_isozymes)),
        rid -> Symbol.(keys(reaction_isozymes[string(rid)])),
    )

    # allocate variables for everything (nb. += wouldn't associate right here)
    constraints =
        constraints +
        :fluxes_forward^unsigned_positive_contribution_variables(fluxes) +
        :fluxes_reverse^unsigned_negative_contribution_variables(fluxes) +
        :isozyme_forward_amounts^isozyme_amounts +
        :isozyme_reverse_amounts^isozyme_amounts +
        :gene_product_amounts^C.variables(keys = gene_ids, bounds = C.Between(0, Inf))

    # connect all parts with constraints
    constraints =
        constraints *
        :directional_flux_balance^sign_split_constraints(
            positive = constraints.fluxes_forward,
            negative = constraints.fluxes_reverse,
            signed = fluxes,
        ) *
        :isozyme_flux_forward_balance^isozyme_flux_constraints(
            constraints.isozyme_forward_amounts,
            constraints.fluxes_forward,
            (rid, isozyme) -> maybemap(
                x -> x.kcat_forward,
                maybeget(reaction_isozymes, string(rid), string(isozyme)),
            ),
        ) *
        :isozyme_flux_reverse_balance^isozyme_flux_constraints(
            constraints.isozyme_reverse_amounts,
            constraints.fluxes_reverse,
            (rid, isozyme) -> maybemap(
                x -> x.kcat_reverse,
                maybeget(reaction_isozymes, string(rid), string(isozyme)),
            ),
        ) *
        :gene_product_isozyme_balance^gene_product_isozyme_constraints(
            constraints.gene_product_amounts,
            (constraints.isozyme_forward_amounts, constraints.isozyme_reverse_amounts),
            (rid, isozyme) -> maybemap(
                x -> [(Symbol(k), v) for (k, v) in x.gene_product_stoichiometry],
                maybeget(reaction_isozymes, string(rid), string(isozyme)),
            ),
        ) *
        :gene_product_capacity^(
            capacity isa Float64 ?
            C.Constraint(
                value = sum(
                    gpa.value * gene_product_molar_masses[String(gp)] for
                    (gp, gpa) in constraints.gene_product_amounts
                ),
                bound = C.Between(0, capacity),
            ) :
            C.ConstraintTree(
                Symbol(id) => C.Constraint(
                    value = sum(
                        constraints.gene_product_amounts[Symbol(gp)].value *
                        gene_product_molar_masses[gp] for gp in gps
                    ),
                    bound = C.Between(0, limit),
                ) for (id, gps, limit) in capacity_limits
            )
        )

    constraints
end

export add_enzyme_constraints
