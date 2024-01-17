
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
                            equal_value_constraint(i.value * stoi, gene_product_amounts[gp])
                    end
                end
            end
        end
    end
    res
end

export gene_product_isozyme_constraints
