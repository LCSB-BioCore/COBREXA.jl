
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

"""
$(TYPEDSIGNATURES)

Summarize a dictionary of flux dictionaries obtained eg. from
[`flux_variability_analysis_dict`](@ref). The simplified summary representation
is useful for pretty-printing and easily showing the most important results.

Internally this function uses [`looks_like_biomass_reaction`](@ref) and
[`looks_like_exchange_reaction`](@ref). The corresponding keyword arguments are
passed to these functions. Use this if your model has an uncommon naming of
reactions.

# Example
```
julia> sol = flux_variability_analysis_dict(model, Gurobi.Optimizer; bounds = objective_bounds(0.99))
julia> flux_res = flux_variability_summary(sol)
Biomass                     Lower bound   Upper bound
  BIOMASS_Ecoli_core_w_GAM: 0.8652        0.8652
Exchange
  EX_h2o_e:                 28.34         28.34
  EX_co2_e:                 22.0377       22.0377
  EX_o2_e:                  -22.1815      -22.1815
  EX_h_e:                   17.3556       17.3556
  EX_glc__D_e:              -10.0         -10.0
  EX_nh4_e:                 -4.8448       -4.8448
  EX_pi_e:                  -3.2149       -3.2149
  EX_for_e:                 0.0           0.0
  ...                       ...           ...
```
"""
function flux_variability_summary(
    flux_result::Tuple{Dict{String,Dict{String,Float64}},Dict{String,Dict{String,Float64}}};
    exclude_exchanges = false,
    exchange_prefixes = constants.exchange_prefixes,
    biomass_strings = constants.biomass_strings,
    exclude_biomass = false,
)
    isnothing(flux_result) && return FluxVariabilitySummary()

    rxn_ids = keys(flux_result[1])
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

    biomass_fluxes = Dict{String,Vector{Maybe{Float64}}}()
    for rxn_id in bmasses
        lb = isnothing(flux_result[1][rxn_id]) ? nothing : flux_result[1][rxn_id][rxn_id]
        ub = isnothing(flux_result[2][rxn_id]) ? nothing : flux_result[2][rxn_id][rxn_id]
        biomass_fluxes[rxn_id] = [lb, ub]
    end

    ex_rxn_fluxes = Dict{String,Vector{Maybe{Float64}}}()
    for rxn_id in ex_rxns
        lb = isnothing(flux_result[1][rxn_id]) ? nothing : flux_result[1][rxn_id][rxn_id]
        ub = isnothing(flux_result[2][rxn_id]) ? nothing : flux_result[2][rxn_id][rxn_id]
        ex_rxn_fluxes[rxn_id] = [lb, ub]
    end

    return FluxVariabilitySummary(biomass_fluxes, ex_rxn_fluxes)
end
