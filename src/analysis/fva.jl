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

function fluxVariabilityAnalysis(model::LinearModel, reactions::Array{Int64, 1}, optimizer)
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
