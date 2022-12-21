"""
$(TYPEDSIGNATURES)

Extract the solution of a specific `community_member` from `opt_model`, which is
a solved optimization model built from the `community_model`. Removes the
`community_member` prefix in the string ids of the returned dictionary.
"""
get_solution(community_model::BalancedGrowthCommunityModel, opt_model, community_member::CommunityMember) = is_solved(opt_model) ? Dict(
    string(last(split(rid, community_member.id*"#"))) => val for (rid, val) in zip(reactions(community_model), reaction_variables(community_model)' * value.(opt_model[:x])) if startswith(rid, community_member.id*"#")
) : nothing


"""
$(TYPEDSIGNATURES)

Extract the environmental exchanges from `opt_model`, which is a solved
optimization model built from the `community_model`.
"""
get_environmental_exchanges(community_model::BalancedGrowthCommunityModel, opt_model) = is_solved(opt_model) ? Dict(
    rid => val for (rid, val) in zip(reactions(community_model), reaction_variables(community_model)' * value.(opt_model[:x])) if !any(startswith(rid, cm.id*"#") for cm in community_model.members)
) : nothing

