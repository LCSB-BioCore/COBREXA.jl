
"""
function make_optimization_model(
    model::MetabolicModel,
    optimizer;
    sense = MOI.MAX_SENSE,
)

Convert CoreModel to the JuMP model, place objectives and the equality
constraint.
"""
function make_optimization_model(model::MetabolicModel, optimizer; sense = MOI.MAX_SENSE)
    m, n = size(stoichiometry(model))
    xl, xu = bounds(model)

    optimization_model = COBREXA.JuMP.Model(optimizer)
    @variable(optimization_model, x[i = 1:n])
    @objective(optimization_model, sense, objective(model)' * x)
    @constraint(optimization_model, mb, stoichiometry(model) * x .== balance(model)) # mass balance
    @constraint(optimization_model, lbs, xl .<= x) # lower bounds
    @constraint(optimization_model, ubs, x .<= xu) # upper bounds

    return optimization_model
end

"""
    optimize_model(
        model::MetabolicModel,
        optimizer;
        sense = MOI.MIN_SENSE,
    )

Use JuMP to solve an instance of CoreModel
"""
function optimize_model(model::MetabolicModel, optimizer; sense = MOI.MIN_SENSE)
    optimization_model = make_optimization_model(model, optimizer; sense = sense)
    COBREXA.JuMP.optimize!(optimization_model)
    return optimization_model
end
