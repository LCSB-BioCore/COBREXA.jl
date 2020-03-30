import ..CobraLP
using JuMP
using GLPK

function solveLP(model::CobraLP)
   m, n = size(model.S)

   optimization_model = JuMP.Model(GLPK.Optimizer)
   @variable(optimization_model, x[i=1:n], lower_bound=model.lb[i], upper_bound=model.ub[i])
   @objective(optimization_model, Min, model.c' * x)
   @constraint(optimization_model, model.S * x .== model.b)

   JuMP.optimize!(optimization_model)

   return (optimization_model, x)
end
