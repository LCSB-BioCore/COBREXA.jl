"""
$(TYPEDSIGNATURES)

Extract the semantic solution of the variables for a specific `community_member`
from `res`. Removes the `community_member` prefix in the string ids of the
returned dictionary.
"""
function values_community_member_dict(
    semantics::Symbol,
    res::ModelWithResult{<:Model},
    community_member::CommunityMember,
)

    is_solved(res) || return nothing

    val_d = values_community_member_dict(res, community_member)

    (sem_ids, _, sem_vard, _) = Accessors.Internal.semantics(semantics)

    ids = sem_ids(community_member.model)

    Dict(
        id => sum(v * val_d[k] for (k, v) in sem_vard(community_member.model)[id]) for
        id in ids
    )
end

"""
$(TYPEDSIGNATURES)

Extract the solution of the variables for a specific `community_member` from
`res`. Removes the `community_member` prefix in the string ids of the returned
dictionary.
"""
function values_community_member_dict(
    res::ModelWithResult{<:Model},
    community_member::CommunityMember,
)
    is_solved(res) || return nothing
    cm = res.model
    opt_model = res.result

    Dict(
        string(last(split(vid, community_member.id * "#"))) => value(opt_model[:x][k]) for
        (k, vid) in enumerate(variables(cm)) if startswith(vid, community_member.id * "#")
    )
end
