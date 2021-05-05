"""
    flux_variability_analysis(
        model::MetabolicModel,
        reactions::Vector{Int},
        optimizer;
        modifications = [],
        workers = [myid()],
        bounds = z -> (z,z),
        ret = objective_value,
    )::Matrix{Float64}

Flux variability analysis solves a pair of optimization problems in `model` for
each flux listed in `reactions`:
```
 min,max xᵢ
s.t. S x = b
    xₗ ≤ x ≤ xᵤ
     cᵀx ≥ bounds(Z₀)[1]
     cᵀx ≤ bounds(Z₀)[2]
```
where Z₀:= cᵀx₀ is the objective value of an optimal solution of the associated
FBA problem (see [`flux_balance_analysis`](@ref)).

The `bounds` is a user-supplied function that specifies the objective bounds
for the variability optimizations, by default it restricts the flux objective
value to the precise optimum reached in FBA. It can return `-Inf` and `Inf` in
first and second pair to remove the limit. Use [`gamma_bounds`](@ref) and
[`objective_bounds`](@ref) for simple bounds.

`optimizer` must be set to a `JuMP`-compatible optimizer. The computation of
the individual optimization problems is transparently distributed to `workers`
(see `Distributed.workers()`).

`ret` is a function used to extract results from optimized JuMP models of the
individual reactions. More detailed information can be extracted e.g. by
setting it to `m -> (JuMP.objective_value(m), JuMP.value.(m[:x]))`.

Returns a matrix of extracted `ret` values for minima and maxima, of total size
`length(reactions)`×2. The optimizer result status is not checked by default,
instead `ret` function can access the `JuMP.termination_status` of the model
and react accordingly, depending on user decision.
"""
function flux_variability_analysis(
    model::MetabolicModel,
    reactions::Vector{Int},
    optimizer;
    modifications = [],
    workers = [myid()],
    bounds = z -> (z, z),
    ret = objective_value,
)
    if any(reactions .< 1) || any(reactions .> n_reactions(model))
        throw(DomainError(reactions, "Index exceeds number of reactions."))
    end

    Z = bounds(
        objective_value(
            flux_balance_analysis(model, optimizer; modifications = modifications),
        ),
    )

    # store a JuMP optimization model at all workers
    save_model = :(
        begin
            optmodel = $COBREXA.make_optimization_model($model, $optimizer)
            $COBREXA._FVA_add_constraint(optmodel, $(objective(model)), optmodel[:x], $Z)
            optmodel
        end
    )
    map(fetch, save_at.(workers, :cobrexa_parfva_model, Ref(save_model)))
    save_model = nothing # this has some volume, free it again

    # schedule FVA parts parallely using pmap
    fluxes = dpmap(
        rid -> :($COBREXA._FVA_optimize_reaction(cobrexa_parfva_model, $rid, $ret)),
        CachingPool(workers),
        [-reactions reactions],
    )

    # free the data on workers
    map(fetch, remove_from.(workers, :cobrexa_parfva_model))

    return fluxes
end

"""
    flux_variability_analysis(
        model::MetabolicModel,
        optimizer;
        kwargs...
    )

A simpler version of [`flux_variability_analysis`](@ref) that maximizes and minimizes all reactions in the model. Arguments are forwarded.
"""
function flux_variability_analysis(model::MetabolicModel, optimizer; kwargs...)
    n = n_reactions(model)
    return flux_variability_analysis(model, collect(1:n), optimizer; kwargs...)
end

"""
    flux_variability_analysis_dict(
        model::MetabolicModel,
        optimizer;
        kwargs...
    )

A variant of [`flux_variability_analysis`](@ref) that returns the individual
maximized and minimized fluxes of all reactions as two dictionaries (of
dictionaries). All keyword arguments except `ret` are passed through.
"""
function flux_variability_analysis_dict(model::MetabolicModel, optimizer; kwargs...)
    vs = flux_variability_analysis(
        model,
        optimizer;
        kwargs...,
        ret = m -> JuMP.value.(m[:x]),
    )
    rxns = reactions(model)

    return (
        Dict(zip(rxns, [Dict(zip(rxns, fluxes)) for fluxes in vs[:, 1]])),
        Dict(zip(rxns, [Dict(zip(rxns, fluxes)) for fluxes in vs[:, 2]])),
    )
end

"""
    gamma_bounds(gamma)

A bounds-generating function for [`flux_variability_analysis`](@ref) that
limits the objective value to be at least `gamma*Z₀`, as usual in COBRA
packages. Use as the `bounds` argument:
```
flux_variability_analysis(model, some_optimizer; bounds = gamma_bounds(0.9))
```
"""
gamma_bounds(gamma) = z -> (gamma * z, Inf)

"""
    (tolerance) = z -> begin

A bounds-generating function for [`flux_variability_analysis`](@ref) that
limits the objective value to a small multiple of Z₀. Use as `bounds` argument,
similarly to [`gamma_bounds`](@ref).
"""
objective_bounds(tolerance) = z -> begin
    vs = (z * tolerance, z / tolerance)
    (minimum(vs), maximum(vs))
end

"""
    _FVA_add_constraint(model, c, x, Z)

Internal helper function for adding constraints to a model. Exists mainly
because for avoiding namespace problems on remote workers.
"""
function _FVA_add_constraint(model, c, x, Z)
    Z[1] > -Inf && @constraint(model, c' * x >= Z[1])
    Z[2] < Inf && @constraint(model, c' * x <= Z[2])
end

"""
    _FVA_get_opt(model, rid)

Internal helper for creating the optimized model on a remote worker, for
avoiding namespace problems.
"""
function _FVA_optimize_reaction(model, rid, ret)
    sense = rid > 0 ? MOI.MAX_SENSE : MOI.MIN_SENSE
    var = all_variables(model)[abs(rid)]

    @objective(model, sense, var)
    optimize!(model)

    if termination_status(model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        return ret(model)
    else
        return nothing
    end
end
