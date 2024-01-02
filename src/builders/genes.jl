
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
"""
knockout_constraints(; fluxes::C.ConstraintTree, knockout_test) = C.ConstraintTree(
    rxn => C.Constraint(C.value(fluxes[rxn]), C.EqualTo(0.0)) for
    rxn in keys(fluxes) if knockout_test(rxn)
)

"""
$(TYPEDSIGNATURES)
"""
gene_knockouts(;
    fluxes::C.ConstraintTree,
    ko_genes::Vector{String},
    model::A.AbstractFBCModel,
) = knockout_constraints(;
    fluxes,
    knockout_test = rxn -> begin
        maybe_avail = A.reaction_gene_products_available(
            model,
            string(rxn),
            g -> !(g in ko_genes), # not available if knocked out
        )
        isnothing(maybe_avail) ? false : !maybe_avail # negate here because of knockout_constraints
    end,
)

#TODO remove the bang from here, there's no side effect
"""
$(TYPEDSIGNATURES)
"""
knockout!(ctmodel::C.ConstraintTree, ko_genes::Vector{String}, model::A.AbstractFBCModel) =
    ctmodel * :gene_knockouts^gene_knockouts(; fluxes = ctmodel.fluxes, ko_genes, model)

"""
$(TYPEDSIGNATURES)

Pipe-able variant.
"""
knockout!(ko_genes::Vector{String}, model::A.AbstractFBCModel) =
    m -> knockout!(m, ko_genes, model)
