"""
$(TYPEDEF)

A standardized structure used to package models that can easily be combined into
a [`CommunityModel`](@ref).

# Fields
$(TYPEDFIELDS)

# Assumptions
1. Exchange reactions, *all* of which are idenitified in
   `exchange_reaction_ids`, have the form: `A[external] ⟷ ∅` where `A` is a
   metabolite. No other exchanges are allowed. This is not checked, but assumed.
2. There is only one biomass reaction in the model.
"""
Base.@kwdef mutable struct CommunityMember
    "Name of model appended to intracellular reactions and metabolites."
    id::String
    "Underlying model."
    model::AbstractMetabolicModel
    "List of all exchange reactions in model."
    exchange_reaction_ids::Vector{String}
    "ID of biomass reaction."
    biomass_reaction_id::String
end

"""
$(TYPEDEF)

A basic structure representing a community model. All `members` are connected
through `environmental_exchange_reactions` if they possess the corresponding
exchange reaction internally (in `member.exchange_reaction_ids`). The
`environmental_exchange_reactions` field is a dictionary with keys corresponding
to the associated metabolite, as well as the reaction lower and upper bounds.

This model structure stitches together individual member models with
environmental exchange reactions, but does not add any objective. Use the
community wrappers for this. The abundance of each member in the community
weights the environmental exchange balances:
```
env_ex_met1 = abundance_1 * ex_met1_member_1 + abundance_2 * ex_met_member_2 + ...
```
Thus, the environmental exchange reactions will have flux units normalized to
total biomass, but the fluxes of each member in the community will be normalized
to its own biomass.

# Fields
$(TYPEDFIELDS)

# Implementation notes
1. All reactions have the `id` of each respective underlying
    [`CommunityMember`](@ref) appended as a prefix with the delimiter `#`.
    Consequently, exchange reactions of the original model will look like
    `species1#EX_...`.
2. All metabolites have the `id` of each respective underlying
    [`CommunityMember`](@ref) appended as a prefix with the delimiter `#`.
3. All genes have the `id` of the respective underlying
    [`CommunityMember`](@ref) appended as a prefix with the delimiter `#`.
4. `environmental_exchange_reactions` is a superset of all exchange reactions
    contained in each underlying member in the community. Only these reactions
    get joined to each underlying model.
"""
Base.@kwdef mutable struct CommunityModel <: AbstractMetabolicModel
    "Models making up the community."
    members::Vector{CommunityMember}
    "Abundances of each community member."
    abundances::Vector{Float64}
    "Environmental exchange reactions ids mapped to their environmental metabolite, and the reaction bounds."
    environmental_exchange_reactions::Dict{String,Vector{String,Float64,Float64}}
end

function Accessors.variables(cm::CommunityModel)
    rxns = [add_community_prefix(m, rid) for m in cm.members for rid in variables(m.model)]
    env_exs = [rid for rid in keys(cm.environmental_exchange_reactions)]
    return [rxns; env_exs]
end

function Accessors.n_variables(cm::CommunityModel)
    num_model_reactions = sum(n_variables(m.model) for m in cm.members)
    num_env_metabolites = length(cm.environmental_exchange_reactions)
    return num_model_reactions + num_env_metabolites
end

function Accessors.metabolites(cm::CommunityModel)
    mets =
        [add_community_prefix(m, mid) for m in cm.members for mid in metabolites(m.model)]
    return [mets; "ENV_" .* [mid for (mid, _, _) in cm.environmental_exchange_reactions]]
end

function Accessors.n_metabolites(cm::CommunityModel)
    num_model_constraints = sum(n_metabolites(m.model) for m in cm.members)
    num_env_metabolites = length(cm.environmental_exchange_reactions)
    return num_model_constraints + num_env_metabolites
end

Accessors.genes(cm::CommunityModel) =
    [add_community_prefix(m, gid) for m in cm.members for gid in genes(m.model)]

Accessors.n_genes(cm::CommunityModel) = sum(n_genes(m.model) for m in cm.members)

Accessors.balance(cm::CommunityModel) = [
    vcat([balance(m.model) for m in cm.members]...)
    spzeros(length(cm.environmental_exchange_reactions))
]

function Accessors.stoichiometry(cm::CommunityModel)

    model_S = blockdiag([stoichiometry(m.model) for m in cm.members]...)
    model_env = spzeros(size(model_S, 1), length(env_met_ids))

    env_mets = [mid for (mid, _, _) in values(cm.environmental_exchange_reactions)]
    env_rxns = keys(cm.environmental_exchange_reactions)

    env_rows = hcat(
        [
            env_ex_matrix(m, env_mets, env_rxns) .* a for
            (m, a) in zip(cm.members, cm.abundances)
        ]...,
    )

    return [
        model_S model_env
        env_rows -I
    ]
end

function Accessors.bounds(cm::CommunityModel)
    models_lbs = vcat([first(bounds(m.model)) for m in cm.members]...)
    models_ubs = vcat([last(bounds(m.model)) for m in cm.members]...)

    env_lbs = [lb for (_, lb, _) in values(cm.environmental_exchange_reactions)]
    env_ubs = [ub for (_, _, ub) in values(cm.environmental_exchange_reactions)]

    return ([models_lbs; env_lbs], [models_ubs; env_ubs])
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
    nr = length(cm.environmental_exchange_reactions)
    blockdiag(rfs, spdiagm(fill(1, nr)))
end

Accessors.reactions(cm::CommunityModel) = [
    vcat([add_community_prefix.(Ref(m), reactions(m.model)) for m in cm.members]...)
    [rid for rid in keys(cm.environmental_exchange_reactions)]
]

Accessors.n_reactions(cm::CommunityModel) =
    sum(n_reactions(m.model) for m in cm.members) +
    length(cm.environmental_exchange_reactions)

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
