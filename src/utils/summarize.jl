"""
$(TYPEDSIGNATURES)

Print a table with the relevant variables of a solved model. This is useful for
pretty-printing and quickly exploring results.

By default, `semantics = :exchange`. The biomass fluxes are always displayed at
the top level, followed by imported and then exported fluxes. Absolute values of
fluxes smaller than `small_flux_bound` and bigger than `large_flux_bound` are
not printed. Heuristic `biomass_strings` prefixes can be supplied in case the
model uses an esoteric name space cf. [`looks_like_biomass_reaction`](@ref).

In some cases

Returns the full, unaltered `modelwithresult` for further processing.

# Example
```
modelwithresult = flux_balance_analysis(model, GLPK.Optimizer) |> summarize

modelwithresult = flux_balance_analysis(model, GLPK.Optimizer) |> summarize(:reaction)

modelwithresult =
    flux_balance_analysis(model, GLPK.Optimizer) |> summarize(;
        namepace_mapping = (model, rid) ->
            reaction_metabolite_map(model, rid; use_annotation = "bigg.metabolite"),
    )

modelwithresult =
    flux_balance_analysis(model, GLPK.Optimizer) |> summarize(;
        namepace_mapping = (model, rid) ->
            reaction_annotation_map(model, rid; use_annotation = "biocyc"),
    )
```
"""
function summarize(
    modelwithresult::ModelWithResult{<:Model},
    semantics::Symbol = :exchange;
    biomass_strings = constants.biomass_strings,
    small_flux_bound = 1.0 / constants.default_reaction_bound^2,
    large_flux_bound = constants.default_reaction_bound,
    namepace_mapping = nothing,
)
    isnothing(modelwithresult) && return nothing
    ns(x) = isnothing(namepace_mapping) ? x : namepace_mapping(modelwithresult.model, x)

    d = values_dict(semantics, modelwithresult)

    biomass_idxs = findall(
        x ->
            looks_like_biomass_reaction(x; biomass_strings = biomass_strings) ||
                is_sbo_biomass_reaction(modelwithresult.model, x),
        variables(modelwithresult.model),
    )

    biomass_fluxes =
        isnothing(biomass_idxs) ? Float64[] :
        value.(modelwithresult.result[:x][biomass_idxs])
    biomass_ids =
        isnothing(biomass_idxs) ? String[] : variables(modelwithresult.model)[biomass_idxs]

    import_ids = [k for (k, v) in d if -large_flux_bound <= v <= -small_flux_bound]
    import_fluxes = [d[k] for k in import_ids]
    import_idxs = sortperm(import_fluxes)

    export_ids = [k for (k, v) in d if small_flux_bound <= v <= large_flux_bound]
    export_fluxes = [d[k] for k in export_ids]
    export_idxs = sortperm(export_fluxes; rev = true)

    ids = [
        biomass_ids
        ns.(import_ids[import_idxs])
        ns.(export_ids[export_idxs])
    ]
    fluxes = [
        biomass_fluxes
        import_fluxes[import_idxs]
        export_fluxes[export_idxs]
    ]

    data = [
        ids fluxes
    ]

    biomass_high = Highlighter((_, i, _) -> i <= length(biomass_ids), foreground = :yellow)
    import_high = Highlighter(
        (_, i, _) ->
            length(biomass_ids) < i <= length(biomass_ids) + length(import_ids),
        foreground = :blue,
    )
    export_high = Highlighter(
        (_, i, _) ->
            length(biomass_ids) + length(import_ids) <
            i <=
            length(biomass_ids) + length(import_ids) + length(export_ids),
        foreground = :red,
    )

    pretty_table(
        data;
        header = [titlecase(String(semantics)), "Flux"],
        highlighters = (biomass_high, import_high, export_high),
    )

    return modelwithresult # return entire solution
end

summarize(; kwargs...) =
    (modelwithresult::ModelWithResult{<:Model}) ->
        summarize(modelwithresult, :exchange; kwargs...)

summarize(semantics::Symbol; kwargs...) =
    (modelwithresult::ModelWithResult{<:Model}) ->
        summarize(modelwithresult, semantics; kwargs...)
