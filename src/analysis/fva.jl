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
function fluxVariabilityAnalysis(model::LinearModel, optimizer)
   (m, n) = size(model.S)
   return fluxVariabilityAnalysis(model, collect(1:n), optimizer)
end

function fluxVariabilityAnalysis(model::LinearModel, reactions::Vector{Int64}, optimizer)
   (maximum(reactions) > length(model.rxns)) && error("Index exceeds number of reactions.")
   γ = 1.
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

function parFVA_get_minmax(model, rid)
    var = JuMP.all_variables(model)[rid]
    JuMP.@objective(model, MOI.MIN_SENSE, var)
    JuMP.optimize!(model)
    min_flux = JuMP.objective_value(model)

    JuMP.@objective(model, MOI.MAX_SENSE, var)
    JuMP.optimize!(model)
    max_flux = JuMP.objective_value(model)

    [min_flux max_flux]
end

function parFVA(model::LinearModel, reactions::Vector{Int}, optimizer, workers)
    if any(reactions .> length(model.rxns))
        throw(ArgumentError("reactions contain an out-of-bounds index"))
    end

    gamma = 1. #TODO parametrize?
    (optimization_model, x0) = fluxBalanceAnalysis(model::LinearModel, optimizer)
    Z0 = JuMP.objective_value(optimization_model)

    # save the model data to all workers (`Ref` avoids broadcasting over the tuple)
    save_at.(workers, :cobrexa_parfva_data, Ref((model, Z0, optimizer, gamma)))

    # make a JuMP optimization model
    save_at.(workers, :cobrexa_parfva_model, Ref(:(begin
        model, Z0, optimizer, gamma = cobrexa_parfva_data
        optmodel, x = COBREXA.makeOptimizationModel(model, optimizer)
        COBREXA.parFVA_add_constraint(optmodel, model.c, x, Z0, gamma)
        optmodel
    end)))

    # schedule FVA parts parallely using pmap
    fluxes = vcat(dpmap(
        rid -> :(COBREXA.parFVA_get_minmax(cobrexa_parfva_model, $rid)),
        CachingPool(workers), reactions)...)

    # free the data on workers
    fetch.(remove_from.(workers, :cobrexa_parfva_data))
    fetch.(remove_from.(workers, :cobrexa_parfva_model))

    return fluxes
end
