"""
$(TYPEDSIGNATURES)

A helper function to get the exchange metabolites in the order of the listed
exchange reactions of a [`CommunityMember`](@ref).
"""
get_exchange_mets(m::CommunityMember) =
    [first(keys(reaction_stoichiometry(m.model, r))) for r in m.exchange_reaction_ids]

"""
$(TYPEDSIGNATURES)

A helper function to get the unique environmental metabolites.
"""
get_env_mets(cm::BalancedGrowthCommunityModel) =
    unique(hcat([get_exchange_mets(m) for m in cm.members]...))

"""
$(TYPEDSIGNATURES)

A helper function that creates an exchange/environmental variable linking matrix
for community member `m` with `abundance`.
"""
function env_ex_matrix(m::CommunityMember, env_mets)
    mat = spzeros(length(env_mets), size(stoichiometry(m.model), 2))
    for (env_met_idx, env_met) in enumerate(env_mets)
        !(env_met in metabolites(m.model)) && continue
        rex = first(indexin([env_met], get_exchange_mets(m)))
        isnothing(rex) && continue
        ex_ridx = first(indexin([m.exchange_reaction_ids[rex]], reactions(m.model)))
        mat[env_met_idx, ex_ridx] = 1.0
    end
    return mat
end
