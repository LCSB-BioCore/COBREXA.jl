"""
$(TYPEDEF)

A struct used to store summary information about the solution
of a constraint based analysis result.

# Fields
$(TYPEDFIELDS)
"""
struct FluxSummary
    biomass_fluxes::OrderedDict{String,Float64}
    import_fluxes::OrderedDict{String,Float64}
    export_fluxes::OrderedDict{String,Float64}
    unbounded_fluxes::OrderedDict{String,Float64}
end

"""
$(TYPEDSIGNATURES)

A default empty constructor for `FluxSummary`.
"""
function FluxSummary()
    FluxSummary(
        OrderedDict{String,Float64}(),
        OrderedDict{String,Float64}(),
        OrderedDict{String,Float64}(),
        OrderedDict{String,Float64}(),
    )
end

