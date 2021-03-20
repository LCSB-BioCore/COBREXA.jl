"""
Convert LinearModel to the JuMP model, place objectives and the equality
constraint.
"""
function makeOptimizationModel(
    model::LM,
    optimizer;
    sense = MOI.MAX_SENSE,
) where {LM<:AbstractCobraModel}
    m, n = size(stoichiometry(model))
    xl, xu = bounds(model)

    optimization_model = JuMP.Model(optimizer)
    @variable(optimization_model, x[i = 1:n], lower_bound = xl[i], upper_bound = xu[i])
    @objective(optimization_model, sense, objective(model)' * x)
    @constraint(optimization_model, stoichiometry(model) * x .== balance(model))

    return (optimization_model, x)
end

"""
Use JuMP to solve an instance of LinearModel
"""
function solveLP(model::LM, optimizer; sense = MOI.MIN_SENSE) where {LM<:AbstractCobraModel}
    optimization_model, x = makeOptimizationModel(model, optimizer; sense = sense)
    JuMP.optimize!(optimization_model)
    return (optimization_model, x)
end
