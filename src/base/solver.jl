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
    @variable(optimization_model, x[i = 1:n])
    @objective(optimization_model, sense, objective(model)' * x)
    mb = @constraint(optimization_model, mb, stoichiometry(model) * x .== balance(model)) # mass balance
    lbs = @constraint(optimization_model, lbs, xl .<= x) # lower bounds
    ubs = @constraint(optimization_model, ubs, x .<= xu) # upper bounds

    return optimization_model, x, mb, lbs, ubs
end

"""
Use JuMP to solve an instance of LinearModel
"""
function solveLP(model::LM, optimizer; sense = MOI.MIN_SENSE) where {LM<:AbstractCobraModel}
    optimization_model, x, _, _, _ = makeOptimizationModel(model, optimizer; sense = sense)
    JuMP.optimize!(optimization_model)
    return (optimization_model, x)
end
