"""
Flux variability analysis (FVA)

FVA solves the pair of optimization problems for each flux xᵢ
```
min/max xᵢ
s.t. S x = b
     xₗ ≤ x ≤ xᵤ
     cᵀx ≥ γ Z₀
```
where Z₀:= cᵀx₀ is the objective value of an optimal solution to the associated
FBA problem
"""
function fluxVariabilityAnalysis(
    model::LM,
    optimizer;
    gamma::AbstractFloat = 1.0,
) where {LM<:AbstractLinearModel}
    n = nReactions(model)
    return fluxVariabilityAnalysis(model, collect(1:n), optimizer)
end

function fluxVariabilityAnalysis(
    model::LM,
    reactions::Vector{Int},
    optimizer,
    workers = [myid()];
    gamma::AbstractFloat = 1.0,
) where {LM<:AbstractLinearModel}

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
            COBREXA._parFVA_add_constraint(optmodel, $(objective(model)), x, $Z0, $gamma)
            optmodel
        end
    )
    map(fetch, save_at.(workers, :cobrexa_parfva_model, Ref(save_model)))
    save_model = nothing # this has some volume, free it again

    # schedule FVA parts parallely using pmap
    fluxes = dpmap(
        rid -> :(COBREXA._parFVA_get_opt(cobrexa_parfva_model, $rid)),
        CachingPool(workers),
        [-reactions reactions],
    )

    # free the data on workers
    map(fetch, remove_from.(workers, :cobrexa_parfva_model))

    return fluxes
end

function _parFVA_add_constraint(model, c, x, Z0, gamma)
    JuMP.@constraint(model, c' * x ≥ gamma * Z0)
end

function _parFVA_get_opt(model, rid)
    sense = rid > 0 ? MOI.MAX_SENSE : MOI.MIN_SENSE
    var = JuMP.all_variables(model)[abs(rid)]

    JuMP.@objective(model, sense, var)
    JuMP.optimize!(model)
    return JuMP.objective_value(model)
end
