"""
$(TYPEDSIGNATURES)

Shallow copy of a [`CommunityModel`](@ref)
"""
Base.copy(m::CommunityModel) = CommunityModel(;
    members = m.members,
    abundances = m.abundances,
    environmental_links = m.environmental_links,
)

"""
$(TYPEDSIGNATURES)

Shallow copy of a [`CommunityModel`](@ref)
"""
Base.copy(m::CommunityMember) = CommunityMember(;
    id = m.id,
    model = m.id,
    exchange_reaction_ids = m.exchange_reaction_ids,
    biomass_reaction_id = m.biomass_reaction_id,
)

"""
$(TYPEDSIGNATURES)

Shallow copy of a [`EnvironmentalLink`](@ref)
"""
Base.copy(m::EnvironmentalLink) = EnvironmentalLink(;
    reaction_id = m.reaction_id,
    metabolite_id = m.metabolite_id,
    lower_bound = m.lower_bound,
    upper_bound = m.upper_bound,
)

"""
$(TYPEDSIGNATURES)

A helper function that creates an exchange/environmental variable linking matrix
for community member `m`.
"""
function env_ex_matrix(
    m::CommunityMember,
    env_mets::Vector{String},
    env_rxns::Vector{String},
)
    mat = spzeros(length(env_mets), n_variables(m.model))
    for (midx, mid) in enumerate(env_mets)
        mid in metabolites(m.model) || continue # does not have this exchange
        ridx = first(indexin([env_rxns[midx]], variables(m.model))) #  find index of exchange reaction in member
        mat[midx, ridx] = 1.0
    end
    return mat
end

"""
$(TYPEDSIGNATURES)

A helper function to find the index of the appropriate model. Assumes each `id`
is delimited by `#` that separates the model ID prefix and the original id.
"""
function access_community_member(
    cm::CommunityModel,
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
