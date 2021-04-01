"""
Convert LinearModel to the JuMP model, place objectives and the equality
constraint.
"""
function make_optimization_model(
    model::LM,
    optimizer;
    sense = MOI.MAX_SENSE,
) where {LM<:MetabolicModel}
    m, n = size(stoichiometry(model))
    xl, xu = bounds(model)

    optimization_model = JuMP.Model(optimizer)
    @variable(optimization_model, x[i = 1:n])
    @objective(optimization_model, sense, objective(model)' * x)
    @constraint(optimization_model, mb, stoichiometry(model) * x .== balance(model)) # mass balance
    @constraint(optimization_model, lbs, xl .<= x) # lower bounds
    @constraint(optimization_model, ubs, x .<= xu) # upper bounds

    return optimization_model
end

"""
Use JuMP to solve an instance of LinearModel
"""
function optimize_model(
    model::LM,
    optimizer;
    sense = MOI.MIN_SENSE,
) where {LM<:MetabolicModel}
    optimization_model = make_optimization_model(model, optimizer; sense = sense)
    x = optimization_model[:x]
    JuMP.optimize!(optimization_model)
    return optimization_model
end
