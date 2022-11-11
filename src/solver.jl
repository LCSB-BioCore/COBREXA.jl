
"""
    module Solver

Interface of COBREXA to JuMP solvers; mainly recreation of the
`AbstractMetabolicModel`s into JuMP optimization models.

# Exports
$(EXPORTS)
"""
module Solver
using ..ModuleTools
@dse

using ..Types
using ..Accessors
using JuMP

"""
$(TYPEDSIGNATURES)

Convert an [`AbstractMetabolicModel`](@ref) to a JuMP model for  the `optimizer`, but
return it unsolved. The interface used to construct the optimization model occurs
through the use of accessors. The complete problem is:
```
max 0.5 * x' * Q * x + q' * x
s.t. A * x = b
     dₗ ≤ C * x ≤ dᵤ
     xₗ ≤ x ≤ xᵤ
```
where: 
    1. `Q` is returned by [`quadratic_objective`](@ref)
    2. `q` is returned by [`linear_objective`](@ref)
    3. `A` is returned by [`stoichiometry`](@ref)
    4. `b` is returned by [`balance`](@ref)
    5. `C` is returned by [`coupling`](@ref)
    6. `(dₗ, dᵤ)` is returned by [`coupling_bounds`](@ref)
    7. `(xₗ, xᵤ)` is returned by [`bounds`](@ref)

Additional integer constraints can be included only by modifying the returned
optimization model through modifications.
"""
function make_optimization_model(
    model::AbstractMetabolicModel,
    optimizer;
    sense = MAX_SENSE,
)

    precache!(model)
    
    # get model from accessors
    Q = quadratic_objective(model) 
    q = linear_objective(model)
    A = stoichiometry(model)
    b = balance(model)
    C = coupling(model)
    d_lb_ub = coupling_bounds(model)
    x_lb_ub = bounds(model)
    
    _, n = size(A)
    
    opt_model = Model(optimizer)
    
    @variable(opt_model, x[1:n])
    
    if isnothing(Q) && !isnothing(q) # linear only
        @objective(opt_model, sense, q' * x)
    elseif !isnothing(Q) && isnothing(q) # quadratic only
        @objective(opt_model, sense, 0.5 * x' * Q * x )
    elseif !isnothing(Q) && !isnothing(q) # linear and quadratic
        @objective(opt_model, sense, 0.5 * x' * Q * x + q' * x)
    end # can have not objective assigned -> feasibility problem
    
    isnothing(A) || @constraint(opt_model, mass_balance, A * x .== b) 
    
    if !isnothing(C) 
        @constraint(opt_model, coupling_lower, first(d_lb_ub) .<= C * x )
        @constraint(opt_model, coupling_upper, C * x .<= last(d_lb_ub))
    end
    
    if !isnothing(x_lb_ub) 
        @constraint(opt_model, variable_lower, first(x_lb_ub) .<= x)
        @constraint(opt_model, variable_upper, x .<= last(x_lb_ub)) 
    end

    return opt_model
end

"""
$(TYPEDSIGNATURES)

Return `true` if `opt_model` solved successfully (solution is optimal or locally
optimal).  Return `false` if any other termination status is reached.
Termination status is defined in the documentation of `JuMP`.
"""
is_solved(opt_model) = termination_status(opt_model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]

"""
$(TYPEDSIGNATURES)

Shortcut for running JuMP `optimize!` on a model and returning the objective
value, if solved.
"""
function optimize_objective(opt_model)::Maybe{Float64}
    optimize!(opt_model)
    solved_objective_value(opt_model)
end

"""
$(TYPEDSIGNATURES)

Returns vectors of the lower and upper bounds of `opt_model` constraints, where
`opt_model` is a JuMP model constructed by e.g.
[`make_optimization_model`](@ref) or [`flux_balance_analysis`](@ref).
"""
get_optmodel_bounds(opt_model) = (
    [-normalized_rhs(lb) for lb in opt_model[:lbs]],
    [normalized_rhs(ub) for ub in opt_model[:ubs]],
)

"""
$(TYPEDSIGNATURES)

Helper function to set the bounds of a variable in the model. Internally calls
`set_normalized_rhs` from JuMP. If the bounds are set to `nothing`, they will
not be changed.
"""
function set_optmodel_bound!(
    vidx,
    opt_model;
    lower::Maybe{Real} = nothing,
    upper::Maybe{Real} = nothing,
)
    isnothing(lower) || set_normalized_rhs(opt_model[:lbs][vidx], -lower)
    isnothing(upper) || set_normalized_rhs(opt_model[:ubs][vidx], upper)
end

"""
$(TYPEDSIGNATURES)

Returns the current objective value of a model, if solved.

# Example
```
solved_objective_value(flux_balance_analysis(model, ...))
```
"""
solved_objective_value(opt_model)::Maybe{Float64} =
    is_solved(opt_model) ? objective_value(opt_model) : nothing

@export_locals

end
