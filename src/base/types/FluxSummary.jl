"""
struct FluxSummary

A struct used to store information.
"""
struct FluxSummary
    biomass_fluxes :: OrderedDict{String, Float64}
    import_fluxes :: OrderedDict{String, Float64}
    export_fluxes :: OrderedDict{String, Float64}
    unbounded_fluxes :: OrderedDict{String, Float64}
end

"""
    flux_summary(model::MetabolicModel, flux_result::Dict{String, Float64}; 
        exclude_exchanges = false,
        exchange_prefixes = _constants.exchange_prefixes,
        biomass_strings = _constants.biomass_strings,
        exclude_biomass = false,
        small_flux_bound = 1.0/_constants.default_reaction_bound^2,
        large_flux_bound = _constants.default_reaction_bound,
        round_digits = 3,
        keep_unbounded = false,
        )::FluxSummary

Return a `FluxSummary` struct based on the `flux_result` of a constraint based
analysis simulation of `model`. Internally this function uses
[`looks_like_biomass_reaction`](@ref) and
[`looks_like_exchange_reaction`](@ref). The corresponding keyword arguments can
be set in `summarize` and are then passed to these functions if your model has
nonstandard ids for reactions. Fluxes smaller than `small_flux_bound` are not
displayed, while fluxes larger than `large_flux_bound` are only displayed if
`display_unbounded` is `true`. `round_digits` is used to round the displayed
value of the flux.

This function is most useful as a way to generate a struct that has nice pretty
printing of flux results. The resultant struct can also be used in downstream
applications if necessary.

# Example
```
```
"""
function flux_summary(model::MetabolicModel, flux_result::Dict{String, Float64}; 
    exclude_exchanges = false,
    exchange_prefixes = _constants.exchange_prefixes,
    biomass_strings = _constants.biomass_strings,
    exclude_biomass = false,
    small_flux_bound = 1.0/_constants.default_reaction_bound^2,
    large_flux_bound = _constants.default_reaction_bound,
    keep_unbounded = false,
    )
    
    ex_rxns = filter(x -> looks_like_exchange_reaction(x, exclude_biomass=exclude_biomass, biomass_strings=biomass_strings, exchange_prefixes=exchange_prefixes), reactions(model))
    bmasses = filter(x -> looks_like_biomass_reaction(x; exclude_exchanges=exclude_exchanges, exchange_prefixes=exchange_prefixes, biomass_strings=biomass_strings), reactions(model))

    ex_fluxes = [flux_result[k] for k in ex_rxns]
    bmass_fluxes = [flux_result[k] for k in bmasses]

    idx_srt_fluxes = sortperm(ex_fluxes)
    lower_unbounded = [idx for idx in idx_srt_fluxes if ex_fluxes[idx] <= -large_flux_bound]
    upper_unbounded = [idx for idx in idx_srt_fluxes if ex_fluxes[idx] >= large_flux_bound]
    import_fluxes = [idx for idx in idx_srt_fluxes if -large_flux_bound < ex_fluxes[idx] <= -small_flux_bound]
    export_fluxes =  [idx for idx in idx_srt_fluxes if small_flux_bound < ex_fluxes[idx] <= large_flux_bound]

    if keep_unbounded
        return FluxSummary(Dict(k => v for (k, v) in zip(bmasses, bmass_fluxes)), 
        OrderedDict(k => v for (k, v) in zip(ex_rxns[import_fluxes], ex_fluxes[import_fluxes])),
        OrderedDict(k => v for (k, v) in zip(ex_rxns[export_fluxes], ex_fluxes[export_fluxes])),
        OrderedDict(k => v for (k, v) in zip([ex_rxns[lower_unbounded]; ex_rxns[upper_unbounded]], [ex_fluxes[lower_unbounded]; ex_fluxes[upper_unbounded]])),
        )
    else
        return FluxSummary(OrderedDict(k => v for (k, v) in zip(bmasses, bmass_fluxes)), 
        OrderedDict(k => v for (k, v) in zip(ex_rxns[import_fluxes], ex_fluxes[import_fluxes])),
        OrderedDict(k => v for (k, v) in zip(ex_rxns[export_fluxes], ex_fluxes[export_fluxes])),
        OrderedDict{String, Float64}(),
        )
    end
end
