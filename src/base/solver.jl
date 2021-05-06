
"""
    make_optimization_model(
        model::MetabolicModel,
        optimizer;
        sense = MOI.MAX_SENSE,
    )

Convert `MetabolicModel`s to a JuMP model, place objectives and the equality
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

    C = coupling(model) # empty if no coupling
    cl, cu = coupling_bounds(model)
    isempty(C) || @constraint(optimization_model, c_lbs, cl .<= coupling(model) * x) # coupling lower bounds
    isempty(C) || @constraint(optimization_model, c_ubs, coupling(model) * x .<= cu) # coupling upper bounds

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


"""
    is_solved(optmodel)

Return `true` if `optmodel` solved successfully (solution is optimal or locally optimal).
Return `false` if any other termination status is reached. Termination status is defined
in the documentation of `JuMP`.
"""
function is_solved(optmodel)
    COBREXA.JuMP.termination_status(optmodel) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED] ? true : false
end
