"""
$(TYPEDEF)

A structure for representing a single reaction in a [`ObjectModel`](@ref).

# Fields
$(TYPEDFIELDS)

Default constructor requires the `id` keyword to be assigned.
"""
Base.@kwdef mutable struct Reaction
    id::String
    name::Maybe{String} = nothing
    metabolites::Dict{String,Float64} = Dict{String,Float64}()
    lower_bound::Float64 = -constants.default_reaction_bound
    upper_bound::Float64 = constants.default_reaction_bound
    gene_associations::Maybe{GeneAssociations} = nothing
    subsystem::Maybe{String} = nothing
    notes::Notes = Notes()
    annotations::Annotations = Annotations()
end

"""
$(TYPEDSIGNATURES)

A convenience constructor for [`Reaction`](@ref) taking as first argument the id
of the reaction. All other kwargs are forwarded to the type constructor.
"""
Reaction(id; kwargs...) = Reaction(; id, kwargs...)

"""
$(TYPEDSIGNATURES)

Convenience constructor for `Reaction` that generates a reaction constrained to
carry flux only in the forward direction relative to the `metabolites`, which is
a dictionary mapping metabolite ids to stoichiometric coefficients. The
`default_bound` is the value taken to mean infinity in the context of constraint
based models, often this is set to a very high flux value like 1000.

See also: [`Reaction`](@ref), [`ReactionBackward`](@ref), [`ReactionBidirectional`](@ref)
"""
ReactionForward(id::String, metabolites; default_bound = constants.default_reaction_bound) =
    Reaction(id; metabolites = metabolites, lower_bound = 0.0, upper_bound = default_bound)

"""
$(TYPEDSIGNATURES)

Convenience constructor for `Reaction` that generates a reaction constrained to
carry flux only in the backward direction relative to the `metabolites`, which is
a dictionary mapping metabolite ids to stoichiometric coefficients. The
`default_bound` is the value taken to mean infinity in the context of constraint
based models, often this is set to a very high flux value like 1000.

See also: [`Reaction`](@ref), [`ReactionForward`](@ref), [`ReactionBidirectional`](@ref)
"""
ReactionBackward(
    id::String,
    metabolites;
    default_bound = constants.default_reaction_bound,
) = Reaction(id; metabolites = metabolites, lower_bound = -default_bound, upper_bound = 0.0)

"""
$(TYPEDSIGNATURES)

Convenience constructor for `Reaction` that generates a reaction constrained to
carry flux in both the forward and backward direction relative to the
`metabolites`, which is a dictionary mapping metabolite ids to stoichiometric
coefficients. The `default_bound` is the value taken to mean infinity in the
context of constraint based models, often this is set to a very high flux value
like 1000.

See also: [`Reaction`](@ref), [`ReactionForward`](@ref), [`ReactionBackward`](@ref)
"""
ReactionBidirectional(
    id::String,
    metabolites;
    default_bound = constants.default_reaction_bound,
) = Reaction(
    id;
    metabolites = metabolites,
    lower_bound = -default_bound,
    upper_bound = default_bound,
)
