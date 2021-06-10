"""
summarize(model::MetabolicModel, flux_result::Dict{String, Float64}; 
    exclude_exchanges = false,
    exchange_prefixes = _constants.exchange_prefixes,
    biomass_strings = _constants.biomass_strings,
    exclude_biomass = false,
    small_flux_bound = 1.0/_constants.default_reaction_bound^2,
    large_flux_bound = _constants.default_reaction_bound,
    round_digits = 3,
    display_unbounded = false,
    )

Print the `flux_result` of a constraint based analysis simulation of `model`.
Internally this function uses [`looks_like_biomass_reaction`](@ref) and 
[`looks_like_exchange_reaction`](@ref). The corresponding keyword arguments
can be set in `summarize` and are then passed to these functions if your model 
has nonstandard ids for reactions. Fluxes smaller than `small_flux_bound` are 
not displayed, while fluxes larger than `large_flux_bound` are only displayed 
if `display_unbounded` is `true`. `round_digits` is used to round the displayed 
value of the flux.

# Example
```
julia> summarize(model, flux_sol; display_unbounded=true, large_flux_bound=20.0)
Biomass:
        BIOMASS_Ecoli_core_w_GAM: 0.874
Import:
        EX_glc__D_e: -10.0
        EX_nh4_e: -4.765
        EX_pi_e: -3.215
Export:
        EX_h_e: 17.531
Unbounded:
        EX_o2_e: -21.799
        EX_co2_e: 22.81
        EX_h2o_e: 29.176
```
"""
function summarize(model::MetabolicModel, flux_result::Dict{String, Float64}; 
    exclude_exchanges = false,
    exchange_prefixes = _constants.exchange_prefixes,
    biomass_strings = _constants.biomass_strings,
    exclude_biomass = false,
    small_flux_bound = 1.0/_constants.default_reaction_bound^2,
    large_flux_bound = _constants.default_reaction_bound,
    round_digits = 3,
    display_unbounded = false,
    )
    
    ex_rxns = filter(x -> looks_like_exchange_reaction(x, exclude_biomass=exclude_biomass, biomass_strings=biomass_strings, exchange_prefixes=exchange_prefixes), reactions(model))
    bmasses = filter(x -> looks_like_biomass_reaction(x; exclude_exchanges=exclude_exchanges, exchange_prefixes=exchange_prefixes, biomass_strings=biomass_strings), reactions(model))

    ex_fluxes = [flux_result[k] for k in ex_rxns]
    bmass_fluxes = [flux_result[k] for k in bmasses]

    c = Base.text_colors # of course I added colors :)
    println(c[:bold]*c[:light_black], "Biomass:", c[:normal]) # display all the biomass fluxes
    for (k, v) in zip(bmasses, bmass_fluxes)
        println("\t", c[:light_magenta], k, ": ", c[:normal], round(v, digits=round_digits), c[:normal])
    end

    idx_srt_fluxes = sortperm(ex_fluxes)
    lower_unbounded = [idx for idx in idx_srt_fluxes if ex_fluxes[idx] <= -large_flux_bound]
    upper_unbounded = [idx for idx in idx_srt_fluxes if ex_fluxes[idx] >= large_flux_bound]
    import_fluxes = [idx for idx in idx_srt_fluxes if -large_flux_bound < ex_fluxes[idx] <= -small_flux_bound]
    export_fluxes =  [idx for idx in idx_srt_fluxes if small_flux_bound < ex_fluxes[idx] <= large_flux_bound]
    
    println(c[:bold]*c[:light_black],"Import:", c[:normal])
    for (k, v) in zip(ex_rxns[import_fluxes], ex_fluxes[import_fluxes])
        println("\t", c[:light_magenta], k, ": ", c[:normal], round(v, digits=round_digits))
    end

    println(c[:bold]*c[:light_black],"Export:", c[:normal])
    for (k, v) in zip(ex_rxns[export_fluxes], ex_fluxes[export_fluxes])
        println("\t", c[:light_magenta], k, ": ", c[:normal], round(v, digits=round_digits))
    end

    if display_unbounded
        println(c[:bold]*c[:light_black],"Unbounded:", c[:normal])
        for (k, v) in zip([ex_rxns[lower_unbounded]; ex_rxns[upper_unbounded]], [ex_fluxes[lower_unbounded]; ex_fluxes[upper_unbounded]])
            println("\t", c[:light_magenta], k, ": ", c[:normal], round(v, digits=round_digits))
        end
    end

    return nothing
end
