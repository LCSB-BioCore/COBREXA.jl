"""
    FluxVariabilitySummary

A struct used to store summary information about the 
solution of a flux variability analysis result.
"""
struct FluxVariabilitySummary
    biomass_fluxes :: Dict{String, Vector{Union{Float64, Nothing}}}
    exchange_fluxes :: Dict{String, Vector{Union{Float64, Nothing}}} 
end

"""
flux_variability_summary(flux_result::Tuple{Dict{String, Dict{String, Float64}}, Dict{String, Dict{String, Float64}}}; 
    exclude_exchanges = false,
    exchange_prefixes = _constants.exchange_prefixes,
    biomass_strings = _constants.biomass_strings,
    exclude_biomass = false,
    )
Return a 
"""
function flux_variability_summary(flux_result::Tuple{Dict{String, Dict{String, Float64}}, Dict{String, Dict{String, Float64}}}; 
    exclude_exchanges = false,
    exchange_prefixes = _constants.exchange_prefixes,
    biomass_strings = _constants.biomass_strings,
    exclude_biomass = false,
    )

    rxn_ids = keys(flux_result[1])
    ex_rxns = filter(x -> looks_like_exchange_reaction(x, exclude_biomass=exclude_biomass, biomass_strings=biomass_strings, exchange_prefixes=exchange_prefixes), rxn_ids)
    bmasses = filter(x -> looks_like_biomass_reaction(x; exclude_exchanges=exclude_exchanges, exchange_prefixes=exchange_prefixes, biomass_strings=biomass_strings), rxn_ids)

    biomass_fluxes = Dict{String, Vector{Union{Float64, Nothing}}}()
    for rxn_id in bmasses
        lb = isnothing(flux_result[1][rxn_id]) ? nothing : flux_result[1][rxn_id][rxn_id]
        ub = isnothing(flux_result[2][rxn_id]) ? nothing : flux_result[2][rxn_id][rxn_id]
        biomass_fluxes[rxn_id] = [lb, ub]
    end

    ex_rxn_fluxes = Dict{String, Vector{Union{Float64, Nothing}}}() 
    for rxn_id in ex_rxns
        lb = isnothing(flux_result[1][rxn_id]) ? nothing : flux_result[1][rxn_id][rxn_id]
        ub = isnothing(flux_result[2][rxn_id]) ? nothing : flux_result[2][rxn_id][rxn_id]
        ex_rxn_fluxes[rxn_id] = [lb, ub]
    end

    return FluxVariabilitySummary(biomass_fluxes, ex_rxn_fluxes)
end
