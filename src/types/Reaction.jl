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
    kcat_forward::Maybe{Float64} = nothing
    kcat_backward::Maybe{Float64} = nothing
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

Convenience constructor for `Reaction`. The reaction equation is specified using
`metabolites`, which is a dictionary mapping metabolite ids to stoichiometric
coefficients. The direcion of the reaction is set through `dir` which can take
`:bidirectional`, `:forward`, and `:backward` as values. Finally, the
`default_bound` is the value taken to mean infinity in the context of constraint
based models, often this is set to a very high flux value like 1000.
"""
function Reaction(
    id::String,
    metabolites,
    dir = :bidirectional;
    default_bound = constants.default_reaction_bound,
)
    if dir == :forward
        lower_bound = 0.0
        upper_bound = default_bound
    elseif dir == :backward
        lower_bound = -default_bound
        upper_bound = 0.0
    elseif dir == :bidirectional
        lower_bound = -default_bound
        upper_bound = default_bound
    else
        throw(DomainError(dir, "unsupported direction"))
    end
    Reaction(
        id;
        metabolites = metabolites,
        lower_bound = lower_bound,
        upper_bound = upper_bound,
    )
end
