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
function env_ex_member_matrix(
    m::CommunityMember,
    env_mets::Vector{String},
    env_rxns::Vector{String},
)
    idxs = [
        (i, j) for
        (i, j) in enumerate(indexin(env_rxns, variables(m.model))) if !isnothing(j)
    ]
    sparse(
        first.(idxs),
        last.(idxs),
        ones(length(idxs)),
        length(env_mets),
        n_variables(m.model),
    )
end

"""
$(TYPEDSIGNATURES)

A helper function that creates the entire exchange/environmental variable
linking matrix for a community model.
"""
function env_ex_matrix(cm)
    env_mets = [envlink.metabolite_id for envlink in cm.environmental_links]
    env_rxns = [envlink.reaction_id for envlink in cm.environmental_links]

    hcat(
        [
            env_ex_member_matrix(m, env_mets, env_rxns) .* a for
            (m, a) in zip(values(cm.members), cm.abundances)
        ]...,
    )
end

"""
$(TYPEDSIGNATURES)

Variant of [`env_ex_matrix`](@ref) that takes an explicit abundance matrix (used
in solver modifications.)
"""
function env_ex_matrix(cm, abundances)
    env_mets = [envlink.metabolite_id for envlink in cm.environmental_links]
    env_rxns = [envlink.reaction_id for envlink in cm.environmental_links]

    hcat(
        [
            env_ex_member_matrix(m, env_mets, env_rxns) .* a for
            (m, a) in zip(values(cm.members), abundances)
        ]...,
    )
end

"""
$(TYPEDSIGNATURES)

A helper function to find the index of the appropriate model. Assumes each `id`
is delimited by `#` that separates the model ID prefix and the original id.
"""
function access_community_member(
    cm::CommunityModel,
    delim_id::String,
    accessor::Function;
    delim = "#",
    default = nothing,
)
    modelid_nameid = string.(split(delim_id, delim))
    length(modelid_nameid) == 1 && return default # accessor default

    modelid = first(modelid_nameid)
    nameid = last(modelid_nameid)

    accessor(cm.members[modelid].model, nameid)
end

"""
$(TYPEDSIGNATURES)

A helper function to build the `names_lookup` dictionary for a
[`CommunityModel`](@ref).
"""
function build_community_name_lookup(
    members::OrderedDict{String,CommunityMember};
    delim = "#",
)
    accessors = [variables, reactions, metabolites, genes]
    Dict(
        id => Dict(
            Symbol(accessor) =>
                Dict(k => id * delim * k for k in accessor(member.model)) for
            accessor in accessors
        ) for (id, member) in members
    )
end
