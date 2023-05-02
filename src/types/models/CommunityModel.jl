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
    "Underlying model."
    model::AbstractMetabolicModel
    "List of all exchange reactions in model."
    exchange_reaction_ids::Vector{String}
    "ID of biomass reaction."
    biomass_reaction_id::String
end

"""
$(TYPEDEF)

A helper struct used to store the environmental linking information in a
[`CommunityModel`](@ref). The `reaction_id` is an environmental exchange
reaction ID, and `metabolite_id` is the associated metabolite ID. The
`lower_bound` and `upper_bound` are bounds on the flux of the reaction.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct EnvironmentalLink
    "Environmental exchange reaction ID."
    reaction_id::String
    "Environmental metabolite ID."
    metabolite_id::String
    "Exchange reaction lower bound."
    lower_bound::Float64
    "Environmental reaction upper bound."
    upper_bound::Float64
end

"""
$(TYPEDEF)

A basic structure representing a community model. All `members` are connected
through `environmental_links`. The `members` is an `OrderedDict` mapping the
member ID to a [`CommunityMember`](@ref). The `environmental_links` is a vector
of [`EnvironmentalLink`](@ref)s. If a member model possesses any exchange
reaction in `environmental_links`, then it is connected to the associated
environmental exchange reaction. Only the reactions in `environmental_links` are
linked, any other boundary reaction is not constrained in the community model.

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
4. `environmental_links` is a superset of all exchange reactions contained in
    each underlying member in the community. Only these reactions get joined to
    each underlying model.
"""
Base.@kwdef mutable struct CommunityModel <: AbstractMetabolicModel
    "Models making up the community (ID => model)."
    members::OrderedDict{String,CommunityMember}
    "Abundances of each community member."
    abundances::Vector{Float64}
    "Environmental exchange to model exchange linking structure."
    environmental_links::Vector{EnvironmentalLink}
    "A lookup table mapping: model IDs => accessor symbol => names => community names."
    name_lookup::Dict{String,Dict{Symbol,Dict{String,String}}} =
        build_community_name_lookup(members)
end

function Accessors.variable_ids(cm::CommunityModel)
    rxns = [
        cm.name_lookup[id][:variables][vid] for (id, m) in cm.members for
        vid in variable_ids(m.model)
    ]
    env_exs = [envlink.reaction_id for envlink in cm.environmental_links]
    return [rxns; env_exs]
end

function Accessors.variable_count(cm::CommunityModel)
    num_model_reactions = sum(variable_count(m.model) for m in values(cm.members))
    num_env_metabolites = length(cm.environmental_links)
    return num_model_reactions + num_env_metabolites
end

function Accessors.metabolites(cm::CommunityModel)
    mets = [
        cm.name_lookup[id][:metabolites][mid] for (id, m) in cm.members for
        mid in metabolites(m.model)
    ]
    return [mets; "ENV_" .* [envlink.metabolite_id for envlink in cm.environmental_links]]
end

function Accessors.n_metabolites(cm::CommunityModel)
    num_model_constraints = sum(n_metabolites(m.model) for m in values(cm.members))
    num_env_metabolites = length(cm.environmental_links)
    return num_model_constraints + num_env_metabolites
end

Accessors.genes(cm::CommunityModel) =
    [cm.name_lookup[id][:genes][gid] for (id, m) in cm.members for gid in genes(m.model)]

Accessors.n_genes(cm::CommunityModel) = sum(n_genes(m.model) for m in values(cm.members))

Accessors.balance(cm::CommunityModel) = [
    vcat([balance(m.model) for m in values(cm.members)]...)
    spzeros(length(cm.environmental_links))
]

function Accessors.stoichiometry(cm::CommunityModel)

    model_S = blockdiag([stoichiometry(m.model) for m in values(cm.members)]...)
    model_env = spzeros(size(model_S, 1), length(cm.environmental_links))

    env_rows = environment_exchange_stoichiometry(cm)
    env_link = spdiagm(sum(env_rows, dims = 2)[:])

    return [
        model_S model_env
        env_rows -env_link
    ]
end

function Accessors.variable_bounds(cm::CommunityModel)
    models_lbs = vcat([first(variable_bounds(m.model)) for m in values(cm.members)]...)
    models_ubs = vcat([last(variable_bounds(m.model)) for m in values(cm.members)]...)

    env_lbs = [envlink.lower_bound for envlink in cm.environmental_links]
    env_ubs = [envlink.upper_bound for envlink in cm.environmental_links]

    return ([models_lbs; env_lbs], [models_ubs; env_ubs])
end

Accessors.objective(cm::CommunityModel) = spzeros(variable_count(cm))

function Accessors.coupling(cm::CommunityModel)
    coups = blockdiag([coupling(m.model) for m in values(cm.members)]...)
    n = variable_count(cm)
    return [coups spzeros(size(coups, 1), n - size(coups, 2))]
end

Accessors.n_coupling_constraints(cm::CommunityModel) =
    sum(n_coupling_constraints(m.model) for m in values(cm.members))

function Accessors.coupling_bounds(cm::CommunityModel)
    lbs = vcat([first(coupling_bounds(m.model)) for m in values(cm.members)]...)
    ubs = vcat([last(coupling_bounds(m.model)) for m in values(cm.members)]...)
    return (lbs, ubs)
end

function Accessors.reaction_variables(model::CommunityModel)
    nlu_r(id, x) = model.name_lookup[id][:reactions][x]
    nlu_v(id, x) = model.name_lookup[id][:variables][x]
    r_v = Dict{String,Dict{String,Float64}}()
    for (id, m) in model.members
        r_v_m = reaction_variables(m.model)
        for (k, v) in r_v_m
            r_v[nlu_r(id, k)] = Dict(nlu_v(id, kk) => vv for (kk, vv) in v)
        end
    end
    r_v
end

Accessors.reaction_ids(cm::CommunityModel) = [
    vcat(
        [
            [cm.name_lookup[id][:reactions][rid] for rid in reaction_ids(m.model)] for
            (id, m) in cm.members
        ]...,
    )
    [envlink.reaction_id for envlink in cm.environmental_links]
]

Accessors.reaction_count(cm::CommunityModel) =
    sum(reaction_count(m.model) for m in values(cm.members)) +
    length(cm.environmental_links)


Accessors.environmental_exchange_variables(model::CommunityModel) = Dict(
    rid => Dict(rid => 1.0) for
    rid in [envlink.reaction_id for envlink in model.environmental_links]
)

Accessors.environmental_exchange_ids(model::CommunityModel) =
    [envlink.reaction_id for envlink in model.environmental_links]

Accessors.environmental_exchange_count(model::CommunityModel) =
    length(model.environmental_links)

function Accessors.enzyme_variables(model::CommunityModel)
    nlu(id, x) = model.name_lookup[id][:genes][x]
    e_v = Dict{String,Dict{String,Float64}}()
    for (id, m) in model.members
        e_v_m = enzyme_variables(m.model)
        for (k, v) in e_v_m
            e_v[nlu(id, k)] = Dict(nlu(id, kk) => vv for (kk, vv) in v)
        end
    end
    e_v
end

Accessors.enzyme_ids(cm::CommunityModel) = [
    cm.name_lookup[id][:genes][gid] for (id, m) in cm.members for gid in enzyme_ids(m.model)
]

Accessors.enzyme_count(cm::CommunityModel) = sum(enzyme_count(m.model) for m in cm.members)

"""
$(TYPEDSIGNATURES)

Get a mapping of enzyme groups to variables. See [`enzyme_variables`](@ref).
"""
function Accessors.enzyme_group_variables(model::CommunityModel)
    nlu(id, x) = model.name_lookup[id][:genes][x]
    e_g_v = Dict{String,Dict{String,Float64}}()
    for (id, m) in model.members
        e_g_v_m = enzyme_group_variables(m.model)
        for (k, v) in e_g_v_m
            e_g_v[id*"#"*k] = Dict(nlu(id, kk) => vv for (kk, vv) in v)
        end
    end
    e_g_v
end

Accessors.enzyme_group_ids(cm::CommunityModel) =
    [id * "#" * k for (id, m) in cm.members for k in enzyme_group_ids(m.model)]

Accessors.enzyme_group_count(cm::CommunityModel) =
    sum(enzyme_group_count(m.model) for m in cm.members)

#=
This loops implements the rest of the accessors through access_community_member.
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
