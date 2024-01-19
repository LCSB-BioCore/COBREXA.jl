
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

Optionally if `interface` is specified, an "interface block" will be created
within the constraint tree for later use as a "module" in creating bigger
models (such as communities) using [`join_module_constraints`](@ref). The
possible parameter values include:
- `nothing` -- default, no interface is created
- `:sbo` -- the interface gets created from model's SBO annotations)
- `:identifier_prefixes` -- the interface is guesstimated from commonly
  occurring adhoc reaction ID prefixes used in contemporary models
- `:boundary` -- the interface is created from all reactions that either only
  consume or only produce metabolites

Output interface name can be set via `interface_name`.

See [`Configuration`](@ref) for fine-tuning the default interface creation.
"""
function flux_balance_constraints(
    model::A.AbstractFBCModel;
    interface::Maybe{Symbol} = nothing,
    interface_name = :interface,
)
    rxn_strings = A.reactions(model)
    rxns = Symbol.(rxn_strings)
    mets = Symbol.(A.metabolites(model))
    lbs, ubs = A.bounds(model)
    stoi = A.stoichiometry(model)
    bal = A.balance(model)
    obj = A.objective(model)

    constraints = C.ConstraintTree(
        :fluxes^C.variables(keys = rxns, bounds = zip(lbs, ubs)) *
        :flux_stoichiometry^C.ConstraintTree(
            met => C.Constraint(
                value = C.LinearValue(SparseArrays.sparse(row)),
                bound = C.EqualTo(b),
            ) for (met, row, b) in zip(mets, eachrow(stoi), bal)
        ) *
        :objective^C.Constraint(C.LinearValue(SparseArrays.sparse(obj))),
    )

    add_interface(sym, flt) =
        any(flt) && (
            constraints *=
                interface_name^sym^C.ConstraintTree(
                    r => constraints.fluxes[r] for r in rxns[flt]
                )
        )
    if interface == :sbo
        sbod(sbos, rid) = any(in(sbos), get(A.reaction_annotation(model, rid), "sbo", []))
        add_interface(:exchanges, sbod.(Ref(configuration.exchange_sbos), rxn_strings))
        add_interface(:biomass, sbod.(Ref(configuration.biomass_sbos), rxn_strings))
        add_interface(
            :atp_maintenance,
            sbod.(Ref(configuration.atp_maintenance_sbos), rxn_strings),
        )
        add_interface(:demand, sbod.(Ref(configuration.demand_sbos), rxn_strings))
    elseif interface == :identifier_prefixes
        prefixed(ps, s) = any(p -> hasprefix(p, s), ps)
        add_interface(:exchanges, prefixed.(Ref(configuration.exchange_id_prefixes), s))
        add_interface(:biomass, prefixed.(Ref(configuration.biomass_id_prefixes), s))
        add_interface(:atp_maintenance, in.(s, Ref(configuration.atp_maintenance_ids)))
    elseif interface == :boundary
        add_interface(:boundary, ((all(s .<= 0) || all(s .>= 0)) for s in eachcol(stoi)))
    end

    return constraints
end

export flux_balance_constraints

"""
$(TYPEDSIGNATURES)

Build log-concentration-stoichiometry constraints for the `model`, as used e.g.
by [`max_min_driving_force_analysis`](@ref).

The output constraint tree contains a log-concentration variable for each
metabolite in subtree `log_concentrations`. Individual reactions' total
reactant log concentrations (i.e., all log concentrations of actual reactants
minus all log concentrations of products) have their own variables in
`reactant_log_concentrations`. The values are connected by
`log_concentration_stoichiometry`.

Function `concentration_bound` may return a bound for the log-concentration of
a given metabolite (compatible with `ConstraintTrees.Bound`), or `nothing`.
"""
function log_concentration_constraints(
    model::A.AbstractFBCModel;
    concentration_bound = _ -> nothing,
)
    rxns = Symbol.(A.reactions(model))
    mets = Symbol.(A.metabolites(model))
    stoi = A.stoichiometry(model)

    constraints =
        :log_concentrations^C.variables(keys = mets, bounds = concentration_bound.(mets)) +
        :reactant_log_concentrations^C.variables(keys = rxns)

    cs = C.ConstraintTree()

    for (midx, ridx, coeff) in zip(SparseArrays.findnz(stoi)...)
        rid = rxns[ridx]
        value = constraints.log_concentrations[mets[midx]].value * coeff
        if haskey(cs, rid)
            cs[rid].value += value
        else
            cs[rid] = C.Constraint(; value)
        end
    end

    return constraints * :log_concentration_stoichiometry^cs
end

export log_concentration_constraints
