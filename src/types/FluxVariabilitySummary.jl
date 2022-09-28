"""
$(TYPEDEF)

Stores summary information about the result of a flux variability analysis.

# Fields
$(TYPEDFIELDS)
"""
struct FluxVariabilitySummary
    biomass_fluxes::Dict{String,Vector{Maybe{Float64}}}
    exchange_fluxes::Dict{String,Vector{Maybe{Float64}}}
end

"""
$(TYPEDSIGNATURES)

A default empty constructor for [`FluxVariabilitySummary`](@ref).
"""
function FluxVariabilitySummary()
    FluxVariabilitySummary(
        Dict{String,Vector{Maybe{Float64}}}(),
        Dict{String,Vector{Maybe{Float64}}}(),
    )
end
