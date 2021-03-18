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
function fluxVariabilityAnalysis(model::LinearModel, optimizer, γ::AbstractFloat=1.0)
    (m, n) = size(model.S)
    return fluxVariabilityAnalysis(model, collect(1:n), optimizer, γ)
end

function fluxVariabilityAnalysis(model::LinearModel, reactions::Vector{Int64}, optimizer, γ::AbstractFloat=1.0)
    (maximum(reactions) > length(model.rxns)) && error("Index exceeds number of reactions.")
    fluxes = zeros(length(reactions), 2)

    (optimization_model, x₀) = fluxBalanceAnalysis(model::LinearModel, optimizer)
    Z₀ = JuMP.objective_value(optimization_model)
    x = all_variables(optimization_model)
    @constraint(optimization_model, model.c' * x ≥ γ * Z₀)

    for i in eachindex(reactions)
        sense = MOI.MIN_SENSE
        @objective(optimization_model, sense, x[reactions[i]])
        JuMP.optimize!(optimization_model)
        fluxes[i, 1] = JuMP.objective_value(optimization_model)

        sense = MOI.MAX_SENSE
        JuMP.set_objective_sense(optimization_model, sense)
        JuMP.optimize!(optimization_model)
        fluxes[i, 2] = JuMP.objective_value(optimization_model)
    end
    return fluxes
end

function parFVA_add_constraint(model, c, x, Z0, gamma)
    JuMP.@constraint(model, c' * x ≥ gamma * Z0)
end

function parFVA_get_opt(model, rid)
    sense = rid > 0 ? MOI.MAX_SENSE : MOI.MIN_SENSE
    var = JuMP.all_variables(model)[abs(rid)]

    JuMP.@objective(model, sense, var)
    JuMP.optimize!(model)
    return JuMP.objective_value(model)
end

function parFVA(model::LinearModel, reactions::Vector{Int}, optimizer, workers, gamma::AbstractFloat=1.0)
    if any(reactions .> length(model.rxns))
        throw(ArgumentError("reactions contain an out-of-bounds index"))
    end

    (optimization_model, x0) = fluxBalanceAnalysis(model::LinearModel, optimizer)
    Z0 = JuMP.objective_value(optimization_model)

    # make a JuMP optimization model
    map(
        fetch,
        save_at.(
            workers,
            :cobrexa_parfva_model,
            Ref(:(
                begin
                    optmodel, x = COBREXA.makeOptimizationModel($model, $optimizer)
                    COBREXA.parFVA_add_constraint(optmodel, $(model.c), x, $Z0, $gamma)
                    optmodel
                end
            )),
        ),
    )

    # schedule FVA parts parallely using pmap
    fluxes = dpmap(
        rid -> :(COBREXA.parFVA_get_opt(cobrexa_parfva_model, $rid)),
        CachingPool(workers),
        [-reactions reactions],
    )

    # free the data on workers
    map(fetch, remove_from.(workers, :cobrexa_parfva_model))

    return fluxes
end
