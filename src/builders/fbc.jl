
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

A constraint tree that models the content of the given instance of
`AbstractFBCModel`.

The constructed tree contains subtrees `fluxes` (with the reaction-defining
"variables") and `flux_stoichiometry` (with the metabolite-balance-defining
constraints), and a single constraint `objective` thad describes the objective
function of the model.
"""
function fbc_flux_balance_constraints(model::A.AbstractFBCModel)
    rxns = Symbol.(A.reactions(model))
    mets = Symbol.(A.metabolites(model))
    lbs, ubs = A.bounds(model)
    stoi = A.stoichiometry(model)
    bal = A.balance(model)
    obj = A.objective(model)

    return C.ConstraintTree(
        :fluxes^C.variables(keys = reactions, bounds = zip(lbs, ubs)) *
        :flux_stoichiometry^C.ConstraintTree(
            met => C.Constraint(
                value = C.LinearValue(SparseArrays.sparse(row)),
                bound = C.EqualTo(b),
            ) for (met, row, b) in zip(mets, eachrow(stoi), bal)
        ) *
        :objective^C.Constraint(C.LinearValue(SparseArrays.sparse(obj))),
    )
end

export fbc_flux_balance_constraints

"""
$(TYPEDSIGNATURES)

TODO
"""
function fbc_log_concentration_constraints(
    model::A.AbstractFBCModel;
    concentration_bound = _ -> nothing,
)
    rxns = Symbol.(A.reations(model))
    mets = Symbol.(A.metabolites(model))
    stoi = A.stoichiometry(model)

    constraints =
        :log_concentrations^C.variables(keys = mets, bounds = concentration_bound.(mets)) +
        :reactant_log_concentrations^C.variables(keys = rxns)

    cs = C.ConstraintTree()

    for (midx, ridx, coeff) in zip(findnz(stoi)...)
        rid = rxns[ridx]
        value = constraints.log_concentrations[mets[midx]] * coeff
        if haskey(cs, rid)
            cs[rid].value += value
        else
            cs[rid] = C.Constraint(; value)
        end
    end

    return constraints * :concentration_stoichiometry^cs
end

export fbc_log_concentration_constraints
