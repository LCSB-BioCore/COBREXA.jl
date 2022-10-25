"""
$(TYPEDSIGNATURES)

Flux variability analysis solves a pair of optimization problems in `model` for
each flux `f` described in `fluxes`:
```
min,max fᵀxᵢ
s.t. S x = b
    xₗ ≤ x ≤ xᵤ
     cᵀx ≥ bounds(Z₀)[1]
     cᵀx ≤ bounds(Z₀)[2]
```
where Z₀:= cᵀx₀ is the objective value of an optimal solution of the associated
FBA problem (see [`flux_balance_analysis`](@ref) for a related analysis, also
for explanation of the `modifications` argument).

The `bounds` is a user-supplied function that specifies the objective bounds
for the variability optimizations, by default it restricts the flux objective
value to the precise optimum reached in FBA. It can return `-Inf` and `Inf` in
first and second pair to remove the limit. Use [`gamma_bounds`](@ref) and
[`objective_bounds`](@ref) for simple bounds.

`optimizer` must be set to a `JuMP`-compatible optimizer. The computation of
the individual optimization problems is transparently distributed to `workers`
(see `Distributed.workers()`).  The value of Z₀ can be optionally supplied in
argument `optimal_objective_value`, which prevents this function from calling
the non-parallelizable FBA. Separating the single-threaded FBA and
multithreaded variability computation can be used to improve resource
allocation efficiency in many common use-cases.

`ret` is a function used to extract results from optimized JuMP models of the
individual fluxes. By default, it calls and returns the value of
`JuMP.objective_value`. More information can be extracted e.g. by setting it to
a function that returns a more elaborate data structure; such as `m ->
(JuMP.objective_value(m), JuMP.value.(m[:x]))`.

Returns a matrix of extracted `ret` values for minima and maxima, of total size
(`size(fluxes,2)`,2). The optimizer result status is checked with
[`is_solved`](@ref); `nothing` is returned if the optimization failed for any
reason.

# Example
```
model = load_model("e_coli_core.json")
flux_variability_analysis(model, [1, 2, 3, 42], GLPK.optimizer)
```
"""
function flux_variability_analysis(
    model::AbstractMetabolicModel,
    fluxes::SparseMat,
    optimizer;
    modifications = [],
    workers = [myid()],
    optimal_objective_value = nothing,
    bounds = z -> (z, Inf),
    ret = objective_value,
)
    if size(fluxes, 1) != n_reactions(model)
        throw(
            DomainError(
                size(fluxes, 1),
                "Flux matrix size is not compatible with model reaction count.",
            ),
        )
    end

    Z = bounds(
        isnothing(optimal_objective_value) ?
        objective_value(
            flux_balance_analysis(model, optimizer; modifications = modifications),
        ) : optimal_objective_value,
    )

    flux_vector = [fluxes[:, i] for i = 1:size(fluxes, 2)]

    return screen_optmodel_modifications(
        model,
        optimizer;
        common_modifications = vcat(
            modifications,
            [
                (model, opt_model) -> begin
                    Z[1] > -Inf && @constraint(
                        opt_model,
                        objective(model)' * opt_model[:x] >= Z[1]
                    )
                    Z[2] < Inf && @constraint(
                        opt_model,
                        objective(model)' * opt_model[:x] <= Z[2]
                    )
                end,
            ],
        ),
        args = tuple.([flux_vector flux_vector], [MIN_SENSE MAX_SENSE]),
        analysis = (_, opt_model, flux, sense) ->
            _max_variability_flux(opt_model, flux, sense, ret),
        workers = workers,
    )
end

"""
$(TYPEDSIGNATURES)

An overload of [`flux_variability_analysis`](@ref) that explores the fluxes specified by integer indexes
"""
function flux_variability_analysis(
    model::AbstractMetabolicModel,
    flux_indexes::Vector{Int},
    optimizer;
    kwargs...,
)
    if any((flux_indexes .< 1) .| (flux_indexes .> n_fluxes(model)))
        throw(DomainError(flux_indexes, "Flux index out of range"))
    end

    flux_variability_analysis(
        model,
        reaction_flux(model)[:, flux_indexes],
        optimizer;
        kwargs...,
    )
end

"""
$(TYPEDSIGNATURES)

A simpler version of [`flux_variability_analysis`](@ref) that maximizes and
minimizes all declared fluxes in the model. Arguments are forwarded.
"""
flux_variability_analysis(model::AbstractMetabolicModel, optimizer; kwargs...) =
    flux_variability_analysis(model, reaction_flux(model), optimizer; kwargs...)

"""
$(TYPEDSIGNATURES)
A variant of [`flux_variability_analysis`](@ref) that returns the individual
maximized and minimized fluxes as two dictionaries (of dictionaries). All
keyword arguments except `ret` are passed through.

# Example
```
mins, maxs = flux_variability_analysis_dict(
    model,
    Tulip.Optimizer;
    bounds = objective_bounds(0.99),
    modifications = [
        change_optimizer_attribute("IPM_IterationsLimit", 500),
        change_constraint("EX_glc__D_e"; lower_bound = -10, upper_bound = -10),
        change_constraint("EX_o2_e"; lower_bound = 0, upper_bound = 0),
    ],
)
```
"""
function flux_variability_analysis_dict(model::AbstractMetabolicModel, optimizer; kwargs...)
    vs = flux_variability_analysis(
        model,
        optimizer;
        kwargs...,
        ret = sol -> flux_vector(model, sol),
    )
    flxs = fluxes(model)
    dicts = zip.(Ref(flxs), vs)

    return (Dict(flxs .=> Dict.(dicts[:, 1])), Dict(flxs .=> Dict.(dicts[:, 2])))
end

"""
$(TYPEDSIGNATURES)

Internal helper for maximizing reactions in optimization model.
"""
function _max_variability_flux(opt_model, flux, sense, ret)
    @objective(opt_model, sense, sum(flux .* opt_model[:x]))
    optimize!(opt_model)

    is_solved(opt_model) ? ret(opt_model) : nothing
end

"""
$(TYPEDSIGNATURES)

A variant for [`flux_variability_analysis`](@ref) that examines actual
reactions (selected by their indexes in `reactions` argument) instead of whole
fluxes. This may be useful for models where the sets of reactions and fluxes
differ.
"""
function reaction_variability_analysis(
    model::AbstractMetabolicModel,
    reaction_indexes::Vector{Int},
    optimizer;
    kwargs...,
)
    if any((reaction_indexes .< 1) .| (reaction_indexes .> n_reactions(model)))
        throw(DomainError(reaction_indexes, "Flux index out of range"))
    end

    flux_variability_analysis(
        model,
        sparse(
            reaction_indexes,
            1:length(reaction_indexes),
            1.0,
            n_reactions(model),
            length(reaction_indexes),
        ),
        optimizer;
        kwargs...,
    )
end

"""
$(TYPEDSIGNATURES)

Shortcut for [`reaction_variability_analysis`](@ref) that examines all reactions.
"""
reaction_variability_analysis(model::AbstractMetabolicModel, optimizer; kwargs...) =
    reaction_variability_analysis(
        model,
        collect(1:n_reactions(model)),
        optimizer;
        kwargs...,
    )
