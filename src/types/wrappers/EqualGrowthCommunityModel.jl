"""
$(TYPEDEF)

A wrapper around [`CommunityModel`](@ref) that returns a community model where
the growth rates of all members are constrained to be equal to
`community_objective_id`, which is the community growth rate. The objective of
the resultant model is set to this `community_objective_id`. 

# Notes
1. No biomass metabolite exists (and none are created).

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct EqualGrowthCommunityModel <: AbstractModelWrapper
    inner::CommunityModel
    community_objective_id::String
end

Accessors.unwrap_model(model::EqualGrowthCommunityModel) = model.inner

Accessors.variables(cm::EqualGrowthCommunityModel) =
    [variables(cm.inner); cm.community_objective_id]

Accessors.n_variables(cm::EqualGrowthCommunityModel) = n_variables(cm.inner) + 1

Accessors.metabolites(cm::EqualGrowthCommunityModel) =
    [metabolites(cm.inner); [m.id for m in cm.inner.members]]

Accessors.n_metabolites(cm::EqualGrowthCommunityModel) =
    n_metabolites(cm.inner) + length(cm.inner.members)

Accessors.balance(cm::EqualGrowthCommunityModel) = [
    balance(cm.inner)
    spzeros(length(cm.inner.members))
]

function Accessors.stoichiometry(cm::EqualGrowthCommunityModel)

    S = stoichiometry(cm.inner)
    obj_col = spzeros(size(S, 1))

    obj_links = spzeros(length(cm.inner.members), size(S, 2))
    for i = 1:length(cm.inner.members)
        m = cm.inner.members[i]
        j = first(
            indexin([add_community_prefix(m, m.biomass_reaction_id)], variables(cm.inner)),
        )
        obj_links[i, j] = 1.0
    end

    obj = -ones(length(cm.inner.members))

    return [
        S obj_col
        obj_links obj
    ]
end

function Accessors.bounds(cm::EqualGrowthCommunityModel)
    lbs, ubs = bounds(cm.inner)
    return ([lbs; 0], [ubs; constants.default_reaction_bound])
end

function Accessors.objective(cm::EqualGrowthCommunityModel)
    vec = spzeros(n_variables(cm)) # overwrite objective
    vec[end] = 1.0
    return vec
end

Accessors.coupling(cm::EqualGrowthCommunityModel) =
    [coupling(cm.inner) spzeros(n_coupling_constraints(cm.inner))]

"""
$(TYPEDSIGNATURES)

Returns a matrix, which when multipled by the solution of a constraints based
problem, yields the semantically meaningful fluxes that correspond to
[`reactions`](@ref).
"""
function Accessors.reaction_variables_matrix(cm::EqualGrowthCommunityModel)
    mtx = reaction_variables_matrix(cm.inner)
    u = spzeros(1, 1)
    u[1] = 1.0
    blockdiag(mtx, u)
end

Accessors.reactions(cm::EqualGrowthCommunityModel) =
    [reactions(cm.inner); cm.community_objective_id]

Accessors.n_reactions(cm::EqualGrowthCommunityModel) = n_reactions(cm.inner) + 1
