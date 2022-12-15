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
        ex_ridx = first(indexin([m.exchange_reaction_ids[rex]], variables(m.model)))
        mat[env_met_idx, ex_ridx] = 1.0
    end
    return mat
end

"""
$(TYPEDSIGNATURES)

A helper function to find the index of the appropriate model. Assumes each `id`
is delimited by `#` that separates the model ID prefix and the original id.
"""
function access_community_member(
    cm::BalancedGrowthCommunityModel,
    id::String,
    accessor::Function;
    default = nothing,
)
    id_split = split(id, "#")
    idx = findfirst(startswith(first(id_split)), m.id for m in cm.members)
    isnothing(idx) && return default # can only access inside community member
    accessor(cm.members[idx].model, string(last(id_split)))
end


"""
$(TYPEDSIGNATURES)

A helper function to add the id of the community member as a prefix to some string.
"""
add_community_prefix(m::CommunityMember, str::String; delim = "#") = m.id * delim * str
