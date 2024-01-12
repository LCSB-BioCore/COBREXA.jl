
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

Create an "interface" block for constraints created by
[`flux_balance_constraints`](@ref) (or any compatible constraint tree).
"""
fbc_boundary_constraints(
    constraints::C.ConstraintTree;
    ignore = _ -> false,
    bound = _ -> nothing,
) = network_boundary_constraints(c.fluxes.c.flux_stoichiometry; ignore, bound)

"""
$(TYPEDSIGNATURES)

Create an "interface" block for a constrained reaction network described by
variables in `fluxes` and flux balance in `balances`. Boundary reactions are
assumed to be the ones that either "only create" or "only remove" the balanced
components (metabolites), i.e, each dot product of the reaction description
with all balance components is either strictly non-negative or non-positive.
"""
function network_boundary_constraints(
    reactions,
    balances;
    ignore = _ -> false,
    bound = _ -> nothing,
)
    #TODO
end

"""
$(TYPEDSIGNATURES)

Linearly scale all bounds in a constraint tree by the `factor`.
"""
function scale_bounds(tree::C.ConstraintTree, factor)
    C.map(tree) do c
        isnothing(c.bound) ? c : C.Constraint(value = c.value, bound = factor * c.bound)
    end
end

"""
$(TYPEDSIGNATURES)

Join multiple modules into a bigger module.

Modules are like usual constraint trees, but contain an explicitly declared
`interface` part, marked properly a tuple (the parameters should form a
dictionary constructor that would generally look such as `:module_name =>
(module, module.interface)`; the second tuple member may also be specified just
by name as e.g. `:interface`, or omitted while relying on `default_interface`).

Interface parts get merged and constrained to create a new interface; networks
are copied intact.
"""
function join_modules(
    ps::Pair...;
    default_in_interface = :interface,
    out_interface = :interface,
    out_balance = :interface_balance,
    ignore = _ -> false,
    bound = _ -> nothing,
)

    prep(id::String, x) = prep(Symbol(id), x)
    prep(id::Symbol, mod::C.ConstraintTree) = prep(id, (mod, default_interface))
    prep(id::Symbol, (mod, interface)::Tuple{C.ConstraintTree,Symbol}) =
        prep(id, mod, mod[interface])
    prep(id::Symbol, (mod, interface)::Tuple{C.ConstraintTree,C.ConstraintTree}) =
        (id^(:network^mod * :interface^interface))

    # collect everything into one huge network (while also renumbering the interfaces)
    modules = sum(prep.(ps))

    # extract a union of all non-ignored interface keys
    interface_keys =
        filter(!ignore, collect(union((keys(m.interface) for (_, m) in modules)...)))

    # remove the interface placeholders and add variables for the new interface
    constraints =
        ConstraintTree(id => (m.network) for (id, m) in modules) +
        out_interface^C.variables(keys = interface_keys, bounds = bound.(interface_keys))

    # join everything with the interrace balance and return
    constraints * C.map(constraints.out_interface) do ic
        C.Constraint(
            value = sum(c.value for mod in modules for (k, c) in mod.interface) - ic.value,
        )
    end
end

"""
$(TYPEDSIGNATURES)

Overload of `join_module_constraints` for general key-value containers.
"""
join_module_constraints(kv) = join_module_constraints(kv...)

# TODO equal_growth must be preserved, should be extended to ratios (using an extra variable)
# TODO same for exchanges (ratio_bounds?)
# TODO probably similar to the distance bounds from MMDF -- unify!

"""
$(TYPEDSIGNATURES)

Helper function to create environmental exchange rections.
"""
function environment_exchange_variables(env_ex_rxns = Dict{String,Tuple{Float64,Float64}}())
    rids = collect(keys(env_ex_rxns))
    lbs_ubs = collect(values(env_ex_rxns))
    C.variables(; keys = Symbol.(rids), bounds = lbs_ubs)
end

export environment_exchange_variables

"""
$(TYPEDSIGNATURES)

Helper function to build a "blank" community model with only environmental exchange reactions.
"""
function build_community_environment(env_ex_rxns = Dict{String,Tuple{Float64,Float64}}())
    C.ConstraintTree(
        :environmental_exchange_reactions => environment_exchange_variables(env_ex_rxns),
    )
end

export build_community_environment

"""
$(TYPEDSIGNATURES)

Helper function to link species specific exchange reactions to the environmental
exchange reactions by weighting them with their abundances.
"""
function link_environmental_exchanges(
    m::C.ConstraintTree,
    member_abundances::Vector{Tuple{Symbol,Float64}};
    on = m.:environmental_exchange_reactions,
    member_fluxes_id = :fluxes,
)
    C.ConstraintTree(
        rid => C.Constraint(
            value = -rxn.value + sum(
                abundance * m[member][member_fluxes_id][rid].value for
                (member, abundance) in member_abundances if
                haskey(m[member][member_fluxes_id], rid);
                init = zero(C.LinearValue),
            ),
            bound = 0.0,
        ) for (rid, rxn) in on
    )
end

export link_environmental_exchanges

"""
$(TYPEDSIGNATURES)

Helper function to set each species growth rate equal to each other.
"""
function equal_growth_rate_constraints(
    member_biomasses::Vector{Tuple{Symbol,C.LinearValue}},
)
    C.ConstraintTree(
        Symbol(bid1, :_, bid2) => C.Constraint(value = bval1 - bval2, bound = 0.0) for
        ((bid1, bval1), (bid2, bval2)) in
        zip(member_biomasses[1:end-1], member_biomasses[2:end])
    )
end

export equal_growth_rate_constraints
