"""
    makeOptimizationModel(
        model::LM,
        optimizer;
        sense = MOI.MAX_SENSE,
    ) where {LM<:AbstractCobraModel}

Convert LinearModel to a JuMP model, and place objectives, flux bounds and
equality "balance" constraint.
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
    solveLP(model::LM, optimizer; sense = MOI.MIN_SENSE) where {LM<:AbstractCobraModel}

Use JuMP to solveLP an instance of LinearModel. Returns a tuple that contains the
new model and a vector of its variables.
"""
function solveLP(model::LM, optimizer; sense = MOI.MIN_SENSE) where {LM<:AbstractCobraModel}
    optimization_model, x = makeOptimizationModel(model, optimizer; sense = sense)
    JuMP.optimize!(optimization_model)
    return (optimization_model, x)
end
