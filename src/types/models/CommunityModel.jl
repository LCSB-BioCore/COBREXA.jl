"""
$(TYPEDEF)

A standardized structure used to package models that can easily be combined into
a [`CommunityModel`](@ref).

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

A basic structure representing a community model. All `members` are connected to
`environmental_exchange_reactions` if they possess the corresponding exchange
reaction internally. This structure stitches together individual models with
environmental exchange reactions, but does not add any objective. Use the
community wrappers for this. The abundance of each member in the community
weights the environmental exchange balance:
```
env_ex_met1 = abundance_1 * ex_met1_member_1 + abundance_2 * ex_met_member_2 + ...
```
Thus, the environmental exchange reactions will have flux units normalized to
total biomass, but the fluxes of each member in the community will be normalized
to its own biomass.

# Fields
$(TYPEDFIELDS)

# Assumptions
1. `environmental_exchange_reactions` is a superset of all exchange reactions
   contained in each underlying member in the community.
2. Each exchanged metabolite in the underlying models, identified through the
   exchange reactions, is associated with an environmental exchange reaction.

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
Base.@kwdef mutable struct CommunityModel <: AbstractMetabolicModel
    "Models making up the community."
    members::Vector{CommunityMember}
    environmental_exchange_reactions::Vector{String}
end

function Accessors.variables(cm::CommunityModel)
    rxns = [add_community_prefix(m, rid) for m in cm.members for rid in variables(m.model)]
    env_exs = ["EX_" * env_met for env_met in get_env_mets(cm)]
    return [rxns; env_exs]
end

function Accessors.n_variables(cm::CommunityModel)
    num_model_reactions = sum(n_variables(m.model) for m in cm.members)
    # assume each env metabolite gets an env exchange
    num_env_metabolites = length(get_env_mets(cm))
    return num_model_reactions + num_env_metabolites
end

function Accessors.metabolites(cm::CommunityModel)
    mets = [add_community_prefix(m, mid) for m in cm.members for mid in metabolites(m.model)]
    return [mets; "ENV_" .* get_env_mets(cm)]
end

function Accessors.n_metabolites(cm::CommunityModel)
    num_model_constraints = sum(n_metabolites(m.model) for m in cm.members)
    # assume each env metabolite gets an env exchange
    num_env_metabolites = length(get_env_mets(cm))
    return num_model_constraints + num_env_metabolites
end

Accessors.genes(cm::CommunityModel) =
    [add_community_prefix(m, gid) for m in cm.members for gid in genes(m.model)]

Accessors.n_genes(cm::CommunityModel) =
    sum(n_genes(m.model) for m in cm.members)

Accessors.balance(cm::CommunityModel) = [
    vcat([balance(m.model) for m in cm.members]...)
    spzeros(length(get_env_mets(cm)))
]

function Accessors.stoichiometry(cm::CommunityModel)
    env_met_ids = get_env_mets(cm)

    model_S = blockdiag([stoichiometry(m.model) for m in cm.members]...)
    model_env = spzeros(size(model_S, 1), length(env_met_ids))
   
    env_rows = hcat([env_ex_matrix(m, env_met_ids) .* m.abundance for m in cm.members]...)
    
    return [
        model_S model_env
        env_rows -I
    ]
end

function Accessors.bounds(cm::CommunityModel)
    models_lbs = vcat([first(bounds(m.model)) for m in cm.members]...)
    models_ubs = vcat([last(bounds(m.model)) for m in cm.members]...)
    
    env_mets = get_env_mets(cm)
    env_lbs = fill(-constants.default_reaction_bound, length(env_mets))
    env_ubs = fill(constants.default_reaction_bound, length(env_mets))
   
    return (
        [models_lbs; env_lbs],
        [models_ubs; env_ubs],
    )
end

Accessors.objective(cm::CommunityModel) = spzeros(n_variables(cm))

function Accessors.coupling(cm::CommunityModel)
    coups = blockdiag([coupling(m.model) for m in cm.members]...)
    n = n_variables(cm)
    return [coups spzeros(size(coups, 1), n - size(coups, 2))]
end

Accessors.n_coupling_constraints(cm::CommunityModel) =
    sum(n_coupling_constraints(m.model) for m in cm.members)

function Accessors.coupling_bounds(cm::CommunityModel)
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
function Accessors.reaction_variables_matrix(cm::CommunityModel)
    # TODO add the non-matrix form!
    rfs = blockdiag([reaction_variables_matrix(m.model) for m in cm.members]...)
    nr = length(get_env_mets(cm))
    blockdiag(rfs, spdiagm(fill(1, nr)))
end

Accessors.reactions(cm::CommunityModel) = [
    vcat([add_community_prefix.(Ref(m), reactions(m.model)) for m in cm.members]...)
    ["EX_" * env_met for env_met in get_env_mets(cm)]
]

Accessors.n_reactions(cm::CommunityModel) =
    sum(n_reactions(m.model) for m in cm.members) + length(get_env_mets(cm))

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
        Accessors.$func(cm::CommunityModel, id::String) =
            access_community_member(cm, id, $func; default = $def)
    end
end
