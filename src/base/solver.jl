
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

    precache!(model)

    m, n = size(stoichiometry(model))
    xl, xu = bounds(model)

    optimization_model = Model(optimizer)
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
    optimize!(optimization_model)
    return optimization_model
end


"""
    is_solved(optmodel)

Return `true` if `optmodel` solved successfully (solution is optimal or locally
optimal).  Return `false` if any other termination status is reached.
Termination status is defined in the documentation of `JuMP`.
"""
function is_solved(optmodel)
    termination_status(optmodel) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED] ? true : false
end

"""
    optimize_objective(optmodel)::Union{Float64,Nothing}

Shortcut for running JuMP `optimize!` on a model and returning the objective
value, if solved.
"""
function optimize_objective(optmodel)::Union{Float64,Nothing}
    optimize!(optmodel)
    if is_solved(optmodel)
        objective_value(optmodel)
    else
        nothing
    end
end

"""
    get_optmodel_bounds(opt_model)

Returns vectors of the lower and upper bounds of `opt_model` constraints, where
`opt_model` is a JuMP model constructed by e.g.
[`make_optimization_model`](@ref) or [`flux_balance_analysis`](@ref).
"""
get_optmodel_bounds(opt_model) = (
    [-normalized_rhs(lb) for lb in opt_model[:lbs]],
    [normalized_rhs(ub) for ub in opt_model[:ubs]],
)

"""
    set_optmodel_bound!(index, optimization_model;
        ub=_constants.default_reaction_rate,
        lb=-_constants.default_reaction_rate)

Helper function to set the bounds of variables.
The JuMP `set_normalized_rhs` function is a little confusing,
so this function simplifies setting constraints. In short, JuMP
uses a normalized right hand side representation of constraints,
which means that lower bounds have their sign flipped. This function
does this for you, so you don't have to remember to do this whenever you
change the constraints.

Just supply the constraint `index` and the JuMP model (`opt_model`) that
will be solved, and the variable's bounds will be set to `ub` and `lb`.
"""
function set_optmodel_bound!(
    vind,
    opt_model;
    ub = _constants.default_reaction_rate,
    lb = -_constants.default_reaction_rate,
)
    set_normalized_rhs(opt_model[:lbs][vind], -lb)
    set_normalized_rhs(opt_model[:ubs][vind], ub)
end
