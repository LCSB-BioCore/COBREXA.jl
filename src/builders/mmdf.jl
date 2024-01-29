
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

Build log-concentration-stoichiometry constraints for the `model`, as used e.g.
by [`max_min_driving_force_analysis`](@ref). If `reaction_subset` is supplied,
then only those reactions, and their associated metabolites, are added to the
resulting model.

The output constraint tree contains a log-concentration variable for each
metabolite in subtree `log_concentrations`. Individual reactions' total reactant
log concentrations (i.e., all log concentrations of actual reactants minus all
log concentrations of products) have their own variables in
`reactant_log_concentrations`. The values are connected by
`log_concentration_stoichiometry`.

Function `concentration_bound` may return a bound for the log-concentration of a
given metabolite (compatible with `ConstraintTrees.Bound`), or `nothing`.
"""
function log_concentration_constraints(
    model::A.AbstractFBCModel;
    reaction_subset = Symbol[],
    concentration_bound = _ -> nothing,
)
    # find only reactions, and metabolites used by reaction_subset if necessary
    rxns = isempty(reaction_subset) ? Symbol.(A.reactions(model)) : reaction_subset
    mets =
        isempty(reaction_subset) ? Symbol.(A.metabolites(model)) :
        unique(
            Symbol(mid) for rid in rxns for
            mid in keys(A.reaction_stoichiometry(model, string(rid)))
        )
    stoi = A.stoichiometry(model)

    constraints =
        :log_concentrations^C.variables(keys = mets, bounds = concentration_bound.(mets)) +
        :reactant_log_concentrations^C.variables(keys = rxns)

    cs = C.ConstraintTree()

    midx_lookup = Dict(indexin(mets, Symbol.(A.metabolites(model))) .=> 1:length(mets)) 
    for (ridx_new, ridx_old) in enumerate(indexin(rxns, Symbol.(A.reactions(model))))
        midxs_old, stoich_coeffs = SparseArrays.findnz(stoi[:, ridx_old])
        rid = rxns[ridx_new]
        for (midx_old, stoich_coeff) in zip(midxs_old, stoich_coeffs)
            midx_new = midx_lookup[midx_old]
            value = constraints.log_concentrations[mets[midx_new]].value * stoich_coeff
            if haskey(cs, rid)
                cs[rid].value += value
            else
                cs[rid] = C.Constraint(; value)
            end
        end
    end

    return constraints * :log_concentration_stoichiometry^cs
end

export log_concentration_constraints
