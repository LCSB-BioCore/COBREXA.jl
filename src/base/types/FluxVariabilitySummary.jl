struct FluxVariabilitySummary
    biomass_fluxes :: Dict{String, Float64}
    import_fluxes :: Dict{String, Float64}
    export_fluxes :: Dict{String, Float64}
    unbounded_fluxes :: Dict{String, Float64}
end

function flux_variability_summary(model::MetabolicModel, flux_result::Dict{String, Dict{String, Float64}}; 
    exclude_exchanges = false,
    exchange_prefixes = _constants.exchange_prefixes,
    biomass_strings = _constants.biomass_strings,
    exclude_biomass = false,
    small_flux_bound = 1.0/_constants.default_reaction_bound^2,
    large_flux_bound = _constants.default_reaction_bound,
    round_digits = 3,
    keep_unbounded = false,
    )

end
