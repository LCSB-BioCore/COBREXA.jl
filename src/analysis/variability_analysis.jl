"""
$(TYPEDSIGNATURES)

Flux variability analysis solves a pair of optimization problems in `model` for
each reaction flux `f` described by [`reactions`](@ref):
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

This is the simplest overload of the function which forwards the arguments to
more complicated ones; see documentation of [`variability_analysis`](@ref) for
details. Arguments `reaction_ids` and `reaction_indexes` may be used to
restrict the analysis to a specified subset of reactions.

# Example
```
model = load_model("e_coli_core.json")
flux_variability_analysis(model, GLPK.optimizer)
```
"""
function flux_variability_analysis(
    model::AbstractMetabolicModel,
    optimizer;
    reaction_ids::Maybe{Vector{String}} = nothing,
    reaction_indexes::Maybe{Vector{Int}} = nothing,
    kwargs...,
)
    variability_analysis(
        :reaction,
        model,
        optimizer;
        ids = reaction_ids,
        indexes = reaction_indexes,
        kwargs...,
    )
end

"""
$(TYPEDSIGNATURES)

A variability analysis over a selected semantics, picking up only objects
specified by IDs or indexes from the selected semantics. For semantics
`:reaction`, this is equivalent to [`flux_variability_analysis`](@ref).
"""
function variability_analysis(
    semantics::Symbol,
    model::AbstractMetabolicModel,
    optimizer;
    ids::Maybe{Vector{String}} = nothing,
    indexes::Maybe{Vector{Int}} = nothing,
    kwargs...,
)
    (sem_ids, n_ids, _, sem_varmtx) = Accessors.Internal.semantics(semantics)

    if isnothing(indexes)
        idxs = if isnothing(ids)
            collect(1:n_ids(model))
        else
            indexin(ids, sem_ids(model))
        end
        any(isnothing.(idxs)) &&
            throw(DomainError(ids[isnothing.(idxs)], "Unknown IDs specified"))
        indexes = Int.(idxs)
    end

    if any((indexes .< 1) .| (indexes .> n_ids(model)))
        throw(DomainError(indexes, "Index out of range"))
    end

    variability_analysis(
        model,
        optimizer;
        directions = sem_varmtx(model)[:, indexes],
        kwargs...,
    )
end

"""
$(TYPEDSIGNATURES)

A generic variability analysis function that maximizes and minimizes the flux
in the directions defined by colums of matrix `directions`.

Argument `modifications` is applied to the model as in
[`flux_balance_analysis`](@ref).

The `bounds` argument is a user-supplied function that specifies the objective
bounds for the variability optimizations as a tuple of minimum and maximum. By
default, it restricts the flux objective value to be greater or equal to the
optimum reached in flux balance analysis. It can return `-Inf` or `Inf` in the
first or second field to completely ignore the limit. Use
[`gamma_bounds`](@ref) and [`objective_bounds`](@ref) for simple bounds.

`optimizer` must be set to a `JuMP`-compatible optimizer. The computation of
the individual optimization problems is transparently distributed to `workers`
(see `Distributed.workers()`).  The value of optimal flux can be optionally
supplied in argument `optimal_objective_value`, which prevents this function
from calling the non-parallelizable flux balance analysis. Separating the
single-threaded FBA and multithreaded variability computation can be used to
improve resource allocation efficiency in many common use-cases.

`ret` is a function used to extract results from optimized JuMP models of the
individual fluxes. By default, it calls and returns the value of
`JuMP.objective_value`. More information can be extracted e.g. by setting it to
a function that returns a more elaborate data structure; such as `m ->
(JuMP.objective_value(m), JuMP.value.(m[:x]))`.

Returns a matrix of extracted `ret` values for minima and maxima, of total size
(`size(directions,2)`,2). The optimizer result status is checked with
[`is_solved`](@ref); `nothing` is returned if the optimization failed for any
reason.
"""
function variability_analysis(
    model::AbstractMetabolicModel,
    optimizer;
    directions::SparseMat = spdiagm(fill(1.0, n_variables(model))),
    modifications = [],
    workers = [myid()],
    optimal_objective_value = nothing,
    bounds = z -> (z, Inf),
    ret = objective_value,
)
    if size(directions, 1) != n_variables(model)
        throw(
            DomainError(
                size(directions, 1),
                "Directions matrix size is not compatible with model variable count.",
            ),
        )
    end

    if isnothing(optimal_objective_value)
        optimal_objective_value =
            flux_balance_analysis(model, optimizer; modifications = modifications) |>
            solved_objective_value
    end
    isnothing(optimal_objective_value) && error("model has no feasible solution for FVA")
    Z = bounds(optimal_objective_value)

    flux_vector = [directions[:, i] for i = 1:size(directions, 2)]

    ModelWithResult(
        model,
        screen_optmodel_modifications(
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
        ),
    )
end

"""
$(TYPEDSIGNATURES)

A variant of [`flux_variability_analysis`](@ref) that returns the individual
maximized and minimized fluxes as two dictionaries of dictionaries. All
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
    # TODO generalize this (requires smart analysis results)
    res = flux_variability_analysis(
        model,
        optimizer;
        kwargs...,
        ret = sol -> values_vec(:reaction, ModelWithResult(model, sol)),
    )
    flxs = reactions(res.model)
    dicts = zip.(Ref(flxs), res.result)

    ModelWithResult(
        res.model,
        (Dict(flxs .=> Dict.(dicts[:, 1])), Dict(flxs .=> Dict.(dicts[:, 2]))),
    )
end

"""
$(TYPEDSIGNATURES)

Internal helper for maximizing reactions in optimization model.
"""
function _max_variability_flux(opt_model, flux, sense, ret)
    @objective(opt_model, sense, sum(flux .* opt_model[:x]))
    optimize!(opt_model)

    # TODO should this get ModelWithResult ?
    is_solved(opt_model) ? ret(opt_model) : nothing
end
