"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)

Default constructor requires the `id` keyword to be assigned.
"""
Base.@kwdef mutable struct Gene
    id::String
    name::Maybe{String} = nothing
    sequence::Maybe{String} = nothing
    product_molar_mass::Maybe{Float64} = nothing
    product_lower_bound::Float64 = 0
    product_upper_bound::Float64 = constants.default_gene_product_bound
    notes::Notes = Notes()
    annotations::Annotations = Annotations()
end

"""
$(TYPEDSIGNATURES)

A convenience constructor for [`Gene`](@ref) taking as first argument the id
of the gene. All other kwargs are forwarded to the type constructor.
"""
Gene(id; kwargs...) = Gene(; id, kwargs...)
