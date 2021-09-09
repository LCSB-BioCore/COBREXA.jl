
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

    optimization_model = Model(optimizer; bridge_constraints = false)
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
function optimize_objective(optmodel)::Maybe{Float64}
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
    set_optmodel_bound!(vidx, opt_model;
        ub::Maybe{Real} = nothing,
        lb::Maybe{Real} = nothing,
    )

Helper function to set the bounds of a variable in the model. Internally calls
`set_normalized_rhs` from JuMP. If the bounds are set to `nothing`, they will
not be changed.
"""
function set_optmodel_bound!(
    vidx,
    opt_model;
    lb::Maybe{Real} = nothing,
    ub::Maybe{Real} = nothing,
)
    isnothing(lb) || set_normalized_rhs(opt_model[:lbs][vidx], -lb)
    isnothing(ub) || set_normalized_rhs(opt_model[:ubs][vidx], ub)
end
