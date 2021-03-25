"""
    fluxVariabilityAnalysis(
        model::LM,
        reactions::Vector{Int},
        optimizer,
        workers = [myid()];
        gamma::AbstractFloat = 1.0,
    )::Matrix{Float64} where {LM<:MetabolicModel}

# Flux variability analysis (FVA)

FVA solves the pair of optimization problems in `model` for each flux xᵢ listed
in `reactions`
```
min/max xᵢ
s.t. S x = b
     xₗ ≤ x ≤ xᵤ
     cᵀx ≥ γ Z₀
```
where Z₀:= cᵀx₀ is the objective value of an optimal solution to the associated
FBA problem.

Internally uses the specified JuMP-compatible `optimizer`, and the work is
scheduled in parallel on `workers`.

Returns a matrix of minima and maxima of size (length(reactions),2).
"""
function fluxVariabilityAnalysis(
    model::LM,
    reactions::Vector{Int},
    optimizer,
    workers = [myid()];
    gamma::AbstractFloat = 1.0,
)::Matrix{Float64} where {LM<:MetabolicModel}

    if any(reactions .< 1) || any(reactions .> nReactions(model))
        throw(DomainError(reactions, "Index exceeds number of reactions."))
    end

    (optimization_model, x0) = fluxBalanceAnalysis(model, optimizer)
    Z0 = JuMP.objective_value(optimization_model)
    optimization_model = nothing # we won't need this one anymore, so free the memory

    # store a JuMP optimization model at all workers
    save_model = :(
        begin
            optmodel, x = COBREXA.makeOptimizationModel($model, $optimizer)
            COBREXA._FVA_add_constraint(optmodel, $(objective(model)), x, $Z0, $gamma)
            optmodel
        end
    )
    map(fetch, save_at.(workers, :cobrexa_parfva_model, Ref(save_model)))
    save_model = nothing # this has some volume, free it again

    # schedule FVA parts parallely using pmap
    fluxes = dpmap(
        rid -> :(COBREXA._FVA_optimize_reaction(cobrexa_parfva_model, $rid)),
        CachingPool(workers),
        [-reactions reactions],
    )

    # free the data on workers
    map(fetch, remove_from.(workers, :cobrexa_parfva_model))

    return fluxes
end

"""
    fluxVariabilityAnalysis(
        model::LM,
        optimizer;
        gamma::AbstractFloat = 1.0,
    ) where {LM<:MetabolicModel}

A simpler version of FVA that maximizes and minimizes all reactions in the model.
"""
function fluxVariabilityAnalysis(
    model::LM,
    optimizer;
    gamma::AbstractFloat = 1.0,
) where {LM<:MetabolicModel}
    n = nReactions(model)
    return fluxVariabilityAnalysis(model, collect(1:n), optimizer; gamma = gamma)
end


"""
    _FVA_add_constraint(model, c, x, Z0, gamma)

Internal helper function for adding constraints to a model. Exists mainly
because for avoiding namespace problems on remote workers.
"""
function _FVA_add_constraint(model, c, x, Z0, gamma)
    JuMP.@constraint(model, c' * x ≥ gamma * Z0)
end

"""
    _FVA_get_opt(model, rid)

Helper for creating the optimized model on a remote worker, for avoiding
namespace problems.
"""
function _FVA_optimize_reaction(model, rid)
    sense = rid > 0 ? MOI.MAX_SENSE : MOI.MIN_SENSE
    var = JuMP.all_variables(model)[abs(rid)]

    JuMP.@objective(model, sense, var)
    JuMP.optimize!(model)
    return JuMP.objective_value(model)
end
