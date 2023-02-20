"""
$(TYPEDEF)

A standardized structure used to package models that can easily be combined into
an [`EqualGrowthCommunityModel`](@ref).

# Fields
$(TYPEDFIELDS)

# Assumptions
1. Exchange reactions, *all* of which are idenitified in `exchange_reaction_ids`,
   have the form: `A[external] ⟷ ∅` where `A` is a metabolite. No other
   exchanges are allowed. This is not checked, but only assumed.
2. There is only one biomass reaction in the model.
"""
Base.@kwdef mutable struct CommunityMember
    "Name of model appended to intracellular reactions and metabolites."
    id::String
    "Abundance of organism."
    abundance::Float64
    "Underlying model."
    model::AbstractMetabolicModel
    "List of all exchange reactions in model."
    exchange_reaction_ids::Vector{String}
    "ID of biomass reaction."
    biomass_reaction_id::String
end

"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)

# Assumptions
1. Each exchanged metabolite in the underlying models, identified through the
   exchange reactions, will get an environmental variable.
2. It is assumed that the same namespace is used to identify unique exchanged
   metabolites.
3. The objective created by this model is the equal growth rate/balanced growth
   objective.

# Implementation notes
1. All reactions have the `id` of each respective underlying
    [`CommunityMember`](@ref) appended as a prefix with the delimiter `#`.
    Consequently, exchange reactions of the original model will look like
    `species1#EX_...`. All exchange environmental reactions have `EX_` as a
    prefix followed by the environmental metabolite id.
2. All metabolites have the `id` of each respective underlying
    [`CommunityMember`](@ref) appended as a prefix with the delimiter `#`. The
    environmental metabolites have no prefix.
3. All genes have the `id` of the respective underlying
    [`CommunityMember`](@ref) appended as a prefix with the delimiter `#`.
"""
Base.@kwdef mutable struct EqualGrowthCommunityModel <: AbstractMetabolicModel
    "Models making up the community."
    members::Vector{CommunityMember}
    "Name of the objective"
    objective_id::String = "equal_growth_rates_biomass_function"
    "Bounds enforced on the environmental exchanged metabolites."
    env_met_flux_bounds::Dict{String,Tuple{Float64,Float64}} =
        Dict{String,Tuple{Float64,Float64}}()
end

function Accessors.variables(cm::EqualGrowthCommunityModel)
    rxns = [add_community_prefix(m, rid) for m in cm.members for rid in variables(m.model)]
    env_exs = ["EX_" * env_met for env_met in get_env_mets(cm)]
    return [rxns; env_exs; cm.objective_id]
end

function Accessors.n_variables(cm::EqualGrowthCommunityModel)
    num_model_reactions = sum(n_variables(m.model) for m in cm.members)
    # assume each env metabolite gets an env exchange
    num_env_metabolites = length(get_env_mets(cm))
    return num_model_reactions + num_env_metabolites + 1 # add 1 for the community biomass
end

function Accessors.metabolites(cm::EqualGrowthCommunityModel)
    mets = [add_community_prefix(m, mid) for m in cm.members for mid in metabolites(m.model)]
    biomass_constraints = [m.id for m in cm.members]
    return [mets; "ENV_" .* get_env_mets(cm); biomass_constraints]
end

function Accessors.n_metabolites(cm::EqualGrowthCommunityModel)
    num_model_constraints = sum(n_metabolites(m.model) for m in cm.members)
    # assume each env metabolite gets an env exchange
    num_env_metabolites = length(get_env_mets(cm))
    return num_model_constraints + num_env_metabolites + length(cm.members)
end

Accessors.genes(cm::EqualGrowthCommunityModel) =
    [add_community_prefix(m, gid) for m in cm.members for gid in genes(m.model)]

Accessors.n_genes(cm::EqualGrowthCommunityModel) =
    sum(n_genes(m.model) for m in cm.members)

Accessors.balance(cm::EqualGrowthCommunityModel) = [
    vcat([balance(m.model) for m in cm.members]...)
    spzeros(length(get_env_mets(cm)))
    spzeros(length(cm.members))
]

function Accessors.stoichiometry(cm::EqualGrowthCommunityModel)
    env_met_ids = get_env_mets(cm)

    model_S = blockdiag([stoichiometry(m.model) for m in cm.members]...)
    model_env = spzeros(size(model_S, 1), length(env_met_ids))
    obj1 = spzeros(size(model_env, 1))

    env_rows = hcat([env_ex_matrix(m, env_met_ids) .* m.abundance for m in cm.members]...)
    obj2 = spzeros(size(env_rows, 1))

    obj_rows = spzeros(length(cm.members), size(model_S, 2))
    for i in 1:length(cm.members)
        m = cm.members[i]
        j = first(indexin([add_community_prefix(m, m.biomass_reaction_id)], variables(cm)))
        obj_rows[i, j] = 1.0
    end
    env_cols = spzeros(length(cm.members), length(env_met_ids))
    obj3 = -ones(length(cm.members))
    
    return [
        model_S model_env obj1
        env_rows -I obj2
        obj_rows env_cols obj3
    ]
end

function Accessors.bounds(cm::EqualGrowthCommunityModel)
    models_lbs = vcat([first(bounds(m.model)) for m in cm.members]...)
    models_ubs = vcat([last(bounds(m.model)) for m in cm.members]...)
    env_lbs = [
        first(get(cm.env_met_flux_bounds, met_id, -constants.default_reaction_bound))
        for met_id in get_env_mets(cm)
    ]
    env_ubs = [
        last(get(cm.env_met_flux_bounds, met_id, constants.default_reaction_bound)) for
        met_id in get_env_mets(cm)
    ]
    return (
        [models_lbs; env_lbs; 0],
        [models_ubs; env_ubs; constants.default_reaction_bound],
    )
end

function Accessors.objective(cm::EqualGrowthCommunityModel)
    vec = spzeros(n_variables(cm))
    vec[end] = 1.0
    return vec
end

function Accessors.coupling(cm::EqualGrowthCommunityModel)
    coups = blockdiag([coupling(m.model) for m in cm.members]...)
    n = n_variables(cm)
    return [coups spzeros(size(coups, 1), n - size(coups, 2))]
end

Accessors.n_coupling_constraints(cm::EqualGrowthCommunityModel) =
    sum(n_coupling_constraints(m.model) for m in cm.members)

function Accessors.coupling_bounds(cm::EqualGrowthCommunityModel)
    lbs = vcat([first(coupling_bounds(m.model)) for m in cm.members]...)
    ubs = vcat([last(coupling_bounds(m.model)) for m in cm.members]...)
    return (lbs, ubs)
end

"""
$(TYPEDSIGNATURES)

Returns a matrix, which when multipled by the solution of a constraints based
problem, yields the semantically meaningful fluxes that correspond to
[`reactions`](@ref).
"""
function Accessors.reaction_variables_matrix(cm::EqualGrowthCommunityModel)
    # TODO add the non-matrix form!
    rfs = blockdiag([reaction_variables_matrix(m.model) for m in cm.members]...)
    nr = length(get_env_mets(cm)) + 1 # env ex + obj
    blockdiag(rfs, spdiagm(fill(1, nr)))
end

Accessors.reactions(cm::EqualGrowthCommunityModel) = [
    vcat([add_community_prefix.(Ref(m), reactions(m.model)) for m in cm.members]...)
    ["EX_" * env_met for env_met in get_env_mets(cm)]
    cm.objective_id
]

Accessors.n_reactions(cm::EqualGrowthCommunityModel) =
    sum(n_reactions(m.model) for m in cm.members) + length(get_env_mets(cm)) + 1

#=
This loops implements the rest of the accssors through access_community_member.
Since most of the environmental reactions are generated programmtically, they
will not have things like annotations etc. For this reason, these methods will
only work if they access something inside the community members.
=#
for (func, def) in (
    (:reaction_gene_associations, nothing),
    (:reaction_subsystem, nothing),
    (:reaction_stoichiometry, nothing),
    (:metabolite_formula, nothing),
    (:metabolite_charge, nothing),
    (:metabolite_compartment, nothing),
    (:reaction_annotations, Dict()),
    (:metabolite_annotations, Dict()),
    (:gene_annotations, Dict()),
    (:reaction_notes, Dict()),
    (:metabolite_notes, Dict()),
    (:gene_notes, Dict()),
    (:reaction_name, nothing),
    (:metabolite_name, nothing),
    (:gene_name, nothing),
)
    @eval begin
        Accessors.$func(cm::EqualGrowthCommunityModel, id::String) =
            access_community_member(cm, id, $func; default = $def)
    end
end
