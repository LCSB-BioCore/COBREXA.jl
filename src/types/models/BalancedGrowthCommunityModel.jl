"""
$(TYPEDEF)

A standardized structure used to package models that can easily be combined into
a [`BalancedGrowthCommunityModel`](@ref).

# Fields
$(TYPEDFIELDS)

# Assumptions
1. Exchange reactions, *all* of which are idenitified in `exchange_reaction_ids`,
   have the form: `A[external] ⟷ ∅` where `A` is a metabolite. No other
   exchanges are allowed.
2. The biomass reaction of a model produces a biomass metabolite called
   `biomass_metabolite_id` with stoichiometry 1.
3. There is only one biomass reaction in the model.
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
    "ID of biomass metabolite."
    biomass_metabolite_id::String
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
   objective. In short, all biomass metabolites are produced at the same rate.
4. The flux units are `mmol X/gDW_total/h` for some metabolite `X`.
"""
Base.@kwdef mutable struct BalancedGrowthCommunityModel <: AbstractMetabolicModel
    "Models making up the community."
    members::Vector{CommunityMember}
    "Name of the objective"
    objective_id::String = "equal_growth_rates_biomass_function"
    "Bounds enforced on the environmental exchanged metabolites."
    env_met_flux_bounds::Dict{String,Tuple{Float64,Float64}} =
        Dict{String,Tuple{Float64,Float64}}()
end

"""
$(TYPEDSIGNATURES)

Return the reactions in `cm`, which is a [`BalancedGrowthCommunityModel`](@ref).
All reactions have the `id` of each respective underlying
[`CommunityMember`](@ref) appended as a prefix with the delimiter `#`.
Consequently, exchange reactions of the original model will look like
`species1#EX_...`. All exchange environmental reactions have `EX_` as a prefix
followed by the environmental metabolite id.
"""
function Accessors.variables(cm::BalancedGrowthCommunityModel)
    rxns = [add_community_prefix(m, rid) for m in cm.members for rid in variables(m.model)]
    env_exs = ["EX_" * env_met for env_met in get_env_mets(cm)]
    return [rxns; env_exs; cm.objective_id]
end

"""
$(TYPEDSIGNATURES)

Return the number of reactions in `cm`, which is a
[`BalancedGrowthCommunityModel`](@ref).
"""
function Accessors.n_variables(cm::BalancedGrowthCommunityModel)
    num_model_reactions = sum(n_variables(m.model) for m in cm.members)
    # assume each env metabolite gets an env exchange
    num_env_metabolites = length(get_env_mets(cm))
    return num_model_reactions + num_env_metabolites + 1 # add 1 for the community biomass
end

"""
$(TYPEDSIGNATURES)

Return the metabolites in `cm`, which is a
[`BalancedGrowthCommunityModel`](@ref). All metabolites have the `id` of each
respective underlying [`CommunityMember`](@ref) appended as a prefix with the
delimiter `#`. The environmental metabolites have no prefix.
"""
function Accessors.constraints(cm::BalancedGrowthCommunityModel)
    mets =
        [add_community_prefix(m, mid) for m in cm.members for mid in metabolites(m.model)]
    return [mets; "ENV_" .* get_env_mets(cm)]
end

"""
$(TYPEDSIGNATURES)

Return the number of metabolites in `cm`, which is a
[`BalancedGrowthCommunityModel`](@ref).
"""
function Accessors.n_constraints(cm::BalancedGrowthCommunityModel)
    num_model_reactions = sum(n_metabolites(m.model) for m in cm.members)
    # assume each env metabolite gets an env exchange
    num_env_metabolites = length(get_env_mets(cm))
    return num_model_reactions + num_env_metabolites
end

"""
$(TYPEDSIGNATURES)

Return the genes in `cm`, which is a [`BalancedGrowthCommunityModel`](@ref). All
genes have the `id` of the respective underlying [`CommunityMember`](@ref)
appended as a prefix with the delimiter `#`.
"""
Accessors.genes(cm::BalancedGrowthCommunityModel) =
    [add_community_prefix(m, gid) for m in cm.members for gid in genes(m.model)]

"""
$(TYPEDSIGNATURES)

Return the number of metabolites in `cm`, which is a [`BalancedGrowthCommunityModel`](@ref).
"""
Accessors.n_genes(cm::BalancedGrowthCommunityModel) =
    sum(n_genes(m.model) for m in cm.members)

"""
$(TYPEDSIGNATURES)

Return the overall stoichiometric matrix for a [`BalancedGrowthCommunityModel`](@ref), built
from the underlying models.
"""
function Accessors.stoichiometry(cm::BalancedGrowthCommunityModel)
    env_mets = get_env_mets(cm)

    model_S = blockdiag([stoichiometry(m.model) for m in cm.members]...)

    zero_rxns = spzeros(size(model_S, 1), length(env_mets))
    obj_rxn = spzeros(size(model_S, 1))
    obj_rxn[indexin(
        [add_community_prefix(m, m.biomass_metabolite_id) for m in cm.members],
        metabolites(cm),
    )] .= [-m.abundance for m in cm.members] # fix units of biomass

    env_met_rows = spzeros(length(env_mets), size(model_S, 2))
    env_met_rows = hcat([env_ex_matrix(m, env_mets) for m in cm.members]...)
    zero_objs = spzeros(length(env_mets))

    return [
        model_S zero_rxns obj_rxn
        env_met_rows -I zero_objs
    ]
end

"""
$(TYPEDSIGNATURES)

Returns the simple variable bounds on `cm`. Assumes the objective can only go
forward at maximum rate `constants.default_reaction_bound`. Note, each bound
from the underlying community members is multiplied by the abundance of that
member.
"""
function Accessors.bounds(cm::BalancedGrowthCommunityModel)
    models_lbs = vcat([first(bounds(m.model)) .* m.abundance for m in cm.members]...)
    models_ubs = vcat([last(bounds(m.model)) .* m.abundance for m in cm.members]...)
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

"""
$(TYPEDSIGNATURES)

Returns the objective of `cm`. This objective is assumed to be the equal growth
rate/balanced growth objective. Consequently, the relation `community_growth *
abundance_species_i = growth_species_i` should hold.
"""
function Accessors.objective(cm::BalancedGrowthCommunityModel)
    vec = spzeros(n_variables(cm))
    vec[end] = 1.0
    return vec
end

"""
$(TYPEDSIGNATURES)

Coupling constraint matrix for a [`BalancedGrowthCommunityModel`](@ref).
"""
function Accessors.coupling(cm::BalancedGrowthCommunityModel)
    coups = blockdiag([coupling(m.model) for m in cm.members]...)
    n = n_variables(cm)
    return [coups spzeros(size(coups, 1), n - size(coups, 2))]
end

"""
$(TYPEDSIGNATURES)

The number of coupling constraints in a [`BalancedGrowthCommunityModel`](@ref).
"""
Accessors.n_coupling_constraints(cm::BalancedGrowthCommunityModel) =
    sum(n_coupling_constraints(m.model) for m in cm.members)

"""
$(TYPEDSIGNATURES)

Coupling bounds for a [`BalancedGrowthCommunityModel`](@ref). Note, each bound
from the underlying community members is multiplied by the abundance of that
member.
"""
function Accessors.coupling_bounds(cm::BalancedGrowthCommunityModel)
    lbs = vcat([first(coupling_bounds(m.model)) .* m.abundance for m in cm.members]...)
    ubs = vcat([last(coupling_bounds(m.model)) .* m.abundance for m in cm.members]...)
    return (lbs, ubs)
end

"""
$(TYPEDSIGNATURES)

Returns a matrix, which when multipled by the solution of a constraints based
problem, yields the semantically meaningful fluxes that correspond to
[`reactions`](@ref).
"""
function Accessors.reaction_variables(cm::BalancedGrowthCommunityModel)
    rfs = blockdiag([reaction_variables(m.model) for m in cm.members]...)
    nr = length(get_env_mets(cm)) + 1 # env ex + obj
    blockdiag(rfs, spdiagm(fill(1, nr)))
end

"""
$(TYPEDSIGNATURES)

Returns the semantically meaningful reactions of the model.
"""
Accessors.reactions(cm::BalancedGrowthCommunityModel) = [
    vcat([add_community_prefix.(Ref(m), reactions(m.model)) for m in cm.members]...)
    ["EX_" * env_met for env_met in get_env_mets(cm)]
    cm.objective_id
]

"""
$(TYPEDSIGNATURES)

Return the semantically meaningful reactions of the model.
"""
Accessors.n_reactions(cm::BalancedGrowthCommunityModel) =
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
    @eval begin # TODO add docstrings somehow
        Accessors.$func(cm::BalancedGrowthCommunityModel, id::String) =
            access_community_member(cm, id, $func; default = $def)
    end
end
