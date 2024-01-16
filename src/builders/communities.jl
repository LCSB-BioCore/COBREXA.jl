
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

Join multiple constraint tree modules with interfaces into a bigger module with
an interface.

Modules are like usual constraint trees, but contain an explicitly declared
`interface` part, marked properly in arguments using e.g. a tuple (the
parameters should form a dictionary constructor that would generally look such
as `:module_name => (module, module.interface)`; the second tuple member may
also be specified just by name as e.g. `:interface`, or omitted while relying
on `default_interface`).

Interface parts get merged and constrained to create a new interface; networks
are intact with disjoint variable sets.

Compatible modules with ready-made interfaces may be created e.g. by
[`flux_balance_constraints`](@ref).
"""
function join_module_constraints(
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

    # first, collect everything into one huge network
    # (while also renumbering the interfaces)
    modules = sum(prep.(ps))

    # fold a union of all non-ignored interface keys
    interface_sum = foldl(modules, init = C.ConstraintTree()) do accs, (_, ms)
        C.imerge(accs, ms) do path, acc, m
            ignore(path) ? missing :
            ismissing(acc) ? C.Constraint(value = m.value) :
            C.Constraint(value = acc.value + m.value)
        end
    end

    # extract the plain networks and add variables for the new interfaces
    constraints =
        ConstraintTree(id => (m.network) for (id, m) in modules) +
        out_interface^C.variables_ifor(bound, interface_sum)

    # join everything with the interrace balance and return
    constraints * C.zip(interface_sum, constraints.out_interface) do sum, out
        C.Constraint(value = sum.value - out.value)
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
