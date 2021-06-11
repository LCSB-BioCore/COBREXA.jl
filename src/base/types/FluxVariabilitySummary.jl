struct FluxVariabilitySummary
    biomass_fluxes :: Dict{String, Vector{Float64}}
    import_fluxes :: Dict{String, Vector{Float64}}
    export_fluxes :: Dict{String, Vector{Float64}}
    unbounded_fluxes :: Dict{String, Vector{Float64}}
end

function flux_variability_summary(model::MetabolicModel, 
    flux_result::Tuple{Dict{String, Dict{String, Float64}}, Dict{String, Dict{String, Float64}}}; 
    exclude_exchanges = false,
    exchange_prefixes = _constants.exchange_prefixes,
    biomass_strings = _constants.biomass_strings,
    exclude_biomass = false,
    small_flux_bound = 1.0/_constants.default_reaction_bound^2,
    large_flux_bound = _constants.default_reaction_bound,
    keep_unbounded = false,
    )

    rxn_ids = keys(flux_result[1])
    ex_rxns = filter(x -> looks_like_exchange_reaction(x, exclude_biomass=exclude_biomass, biomass_strings=biomass_strings, exchange_prefixes=exchange_prefixes), rxn_ids)
    bmasses = filter(x -> looks_like_biomass_reaction(x; exclude_exchanges=exclude_exchanges, exchange_prefixes=exchange_prefixes, biomass_strings=biomass_strings), rxn_ids)


    res = Dict{String, Vector{Float64}}()
    for rxn_id in rxn_ids

        res[rxn_id] = [lb, ub]
    end


end
