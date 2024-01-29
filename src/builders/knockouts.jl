
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

Make a `ConstraintTree` that knocks out fluxes given by the predicate
`knockout_test`. The predicate function is called with a single parameter (the
key of the flux in tree `fluxes`) and must return a boolean. Returning `true`
means that the corresponding flux (usually a reaction flux) will be knocked
out.

Use [`fbc_gene_knockout_constraints`](@ref) to apply gene knockouts easily to
models with `AbstractFBCModel` interface.
"""
knockout_constraints(knockout_test, fluxes::C.ConstraintTree) = C.ConstraintTree(
    id => C.Constraint(C.value(flux), 0) for (id, flux) in fluxes if knockout_test(id)
)

"""
$(TYPEDSIGNATURES)

Make a `ConstraintTree` that simulates a gene knockout of `knockout_genes` in
the `model` and disables corresponding `fluxes` accordingly.

Keys of the fluxes must correspond to the reaction identifiers in the `model`.

`knockout_genes` may be any collection that support element tests using `in`.
Since the test is done many times, a `Set` is a preferred contained for longer
lists of genes.
"""
fbc_gene_knockout_constraints(;
    fluxes::C.ConstraintTree,
    knockout_genes,
    model::A.AbstractFBCModel,
) =
    knockout_constraints(fluxes) do rid
        maybemap(
            !,
            A.reaction_gene_products_available(model, string(rid), !in(knockout_genes)),
        )
    end
