"""
Use JuMP to solve an instance of LinearModel
"""
function solveLP(model::LinearModel, optimizer)
    m, n = size(model.S)

    optimization_model = JuMP.Model(optimizer)
    @variable(
        optimization_model,
        x[i = 1:n],
        lower_bound = model.xl[i],
        upper_bound = model.xu[i]
    )
    @objective(optimization_model, Min, model.c' * x)
    @constraint(optimization_model, model.S * x .== model.b)

    JuMP.optimize!(optimization_model)

    return (optimization_model, x)
end
