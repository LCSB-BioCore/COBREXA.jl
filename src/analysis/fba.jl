"""
Convert LinearModel to the JuMP model, place objectives and the equality
constraint.
"""
function makeOptimizationModel(model::LinearModel, optimizer)
    m,n = size(model.S)

    optimization_model = JuMP.Model(optimizer)
    @variable(optimization_model, x[i=1:n], lower_bound=model.xl[i], upper_bound=model.xu[i])
    @objective(optimization_model, Max, model.c' * x)
    @constraint(optimization_model, model.S * x .== model.b)
    return (optimization_model, x)
end

"""
Flux Balance Analysis
```
max cᵀx
s.t. S x = b
     xₗ ≤ x ≤ xᵤ
```
"""
function fluxBalanceAnalysis(model::LinearModel, optimizer)

   optimization_model, x = makeOptimizationModel(model, optimizer)
   JuMP.optimize!(optimization_model)
   return (optimization_model, x)
   # TODO we might like this to take optimization model and return x, and let
   # the user convert LinearModel to JuMP model. That would make the API much
   # more orthogonal.
end
