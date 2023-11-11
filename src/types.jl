
"""
    Maybe{X}

Type of optional values.
"""
const Maybe{X} = Union{Nothing,X}

"""
$(TYPEDEF)

Information about isozyme composition including subunit stoichiometry and
turnover numbers.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Isozyme
    gene_product_stoichiometry::Dict{String,Float64}
    kcat_forward::Maybe{Float64} = nothing
    kcat_backward::Maybe{Float64} = nothing
end

"""
$(TYPEDSIGNATURES)

A convenience constructor for [`Isozyme`](@ref) that takes a string gene
reaction rule and converts it into the appropriate format. Assumes the
`gene_product_stoichiometry` for each subunit is 1.
"""
Isozyme(gids::Vector{String}; kwargs...) =
    Isozyme(; gene_product_stoichiometry = Dict(gid => 1.0 for gid in gids), kwargs...)
