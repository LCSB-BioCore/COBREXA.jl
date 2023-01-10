"""
$(TYPEDEF)

Information about isozyme composition and annotations justifying the
stoichiometry or turnover numbers.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Isozyme
    stoichiometry::Dict{String,Float64}
    annotation::Annotations = Annotations()
end

"""
$(TYPEDSIGNATURES)

A convenience constructor for [`Isozyme`](@ref) that takes a string gene
reaction rule and converts it into the appropriate format. Assumes the
stoichiometry for each subunit is 1.
"""
Isozyme(gids::Vector{String}; kwargs...) =
    Isozyme(; stoichiometry = Dict(gid => 1.0 for gid in gids), kwargs...)

"""
    const GeneAssociations

An association of genes to a reaction. Each [`Isozyme`](@ref) represents a
distinct enzyme that can catalyze a certain reaction. All the gene products in
an isozyme are required for the enzyme to function. Multiple [`Isozyme`](@ref)s
can catalyze the same reaction, but function independently.
"""
const GeneAssociations = Vector{Isozyme}
