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

function build_community_environment(env_ex_rxns = Dict{String,Tuple{Float64,Float64}}())
    C.ConstraintTree(
        :environmental_exchange_reactions => environment_exchange_variables(env_ex_rxns),
    )
end

export build_community_environment

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

function equal_growth_rate_constraints(member_biomasses::Vector{Tuple{Symbol,C.LinearValue}})
    C.ConstraintTree(
        Symbol(bid1, :_, bid2) => C.Constraint(value = bval1 - bval2, bound = 0.0) for
        ((bid1, bval1), (bid2, bval2)) in
        zip(member_biomasses[1:end-1], member_biomasses[2:end])
    )
end

export equal_growth_rate_constraints
