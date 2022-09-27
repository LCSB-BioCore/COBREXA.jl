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

"""
$(TYPEDSIGNATURES)

Summarize a dictionary of fluxes into small, useful representation of the most
important information contained. Useful for pretty-printing and quickly
exploring the results. Internally this function uses
[`looks_like_biomass_reaction`](@ref) and
[`looks_like_exchange_reaction`](@ref). The corresponding keyword arguments
passed to these functions. Use this if your model has non-standard ids for
reactions. Fluxes smaller than `small_flux_bound` are not stored, while fluxes
larger than `large_flux_bound` are only stored if `keep_unbounded` is `true`.

# Example
```
julia> sol = flux_dict(flux_balance_analysis(model, Tulip.Optimizer))
julia> fr = flux_summary(sol)
Biomass:
  BIOMASS_Ecoli_core_w_GAM: 0.8739
Import:
  EX_o2_e:     -21.7995
  EX_glc__D_e: -10.0
  EX_nh4_e:    -4.7653
  EX_pi_e:     -3.2149
Export:
  EX_h_e:      17.5309
  EX_co2_e:    22.8098
  EX_h2o_e:    29.1758
```
"""
function flux_summary(
    flux_result::Maybe{Dict{String,Float64}};
    exclude_exchanges = false,
    exchange_prefixes = constants.exchange_prefixes,
    biomass_strings = constants.biomass_strings,
    exclude_biomass = false,
    small_flux_bound = 1.0 / constants.default_reaction_bound^2,
    large_flux_bound = constants.default_reaction_bound,
    keep_unbounded = false,
)
    isnothing(flux_result) && return FluxSummary()

    rxn_ids = collect(keys(flux_result))
    ex_rxns = filter(
        x -> looks_like_exchange_reaction(
            x,
            exclude_biomass = exclude_biomass,
            biomass_strings = biomass_strings,
            exchange_prefixes = exchange_prefixes,
        ),
        rxn_ids,
    )
    bmasses = filter(
        x -> looks_like_biomass_reaction(
            x;
            exclude_exchanges = exclude_exchanges,
            exchange_prefixes = exchange_prefixes,
            biomass_strings = biomass_strings,
        ),
        rxn_ids,
    )

    ex_fluxes = [flux_result[k] for k in ex_rxns]
    bmass_fluxes = [flux_result[k] for k in bmasses]

    idx_srt_fluxes = sortperm(ex_fluxes)
    import_fluxes = [
        idx for
        idx in idx_srt_fluxes if -large_flux_bound < ex_fluxes[idx] <= -small_flux_bound
    ]
    export_fluxes = [
        idx for
        idx in idx_srt_fluxes if small_flux_bound < ex_fluxes[idx] <= large_flux_bound
    ]

    if keep_unbounded
        lower_unbounded =
            [idx for idx in idx_srt_fluxes if ex_fluxes[idx] <= -large_flux_bound]
        upper_unbounded =
            [idx for idx in idx_srt_fluxes if ex_fluxes[idx] >= large_flux_bound]
        return FluxSummary(
            OrderedDict(k => v for (k, v) in zip(bmasses, bmass_fluxes)),
            OrderedDict(
                k => v for (k, v) in zip(ex_rxns[import_fluxes], ex_fluxes[import_fluxes])
            ),
            OrderedDict(
                k => v for (k, v) in zip(ex_rxns[export_fluxes], ex_fluxes[export_fluxes])
            ),
            OrderedDict(
                k => v for (k, v) in zip(
                    [ex_rxns[lower_unbounded]; ex_rxns[upper_unbounded]],
                    [ex_fluxes[lower_unbounded]; ex_fluxes[upper_unbounded]],
                )
            ),
        )
    else
        return FluxSummary(
            OrderedDict(k => v for (k, v) in zip(bmasses, bmass_fluxes)),
            OrderedDict(
                k => v for (k, v) in zip(ex_rxns[import_fluxes], ex_fluxes[import_fluxes])
            ),
            OrderedDict(
                k => v for (k, v) in zip(ex_rxns[export_fluxes], ex_fluxes[export_fluxes])
            ),
            OrderedDict{String,Float64}(),
        )
    end
end
