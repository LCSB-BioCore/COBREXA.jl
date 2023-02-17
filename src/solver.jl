
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

Convert `AbstractMetabolicModel`s to a JuMP model, place objectives and the equality
constraint.

Here coupling means inequality constraints coupling multiple variables together.
"""
function make_optimization_model(
    model::AbstractMetabolicModel,
    optimizer;
    sense = MAX_SENSE,
)

    precache!(model)

    m, n = size(stoichiometry(model))
    xl, xu = bounds(model)

    optimization_model = Model(optimizer)
    @variable(optimization_model, x[1:n])
    let obj = objective(model)
        if obj isa AbstractVector
            @objective(optimization_model, sense, obj' * x)
        else
            @objective(optimization_model, sense, x' * obj * [x; 1])
        end
    end
    @constraint(optimization_model, mb, stoichiometry(model) * x .== balance(model)) # mass balance
    @constraint(optimization_model, lbs, xl .<= x) # lower bounds
    @constraint(optimization_model, ubs, x .<= xu) # upper bounds

    C = coupling(model) # empty if no coupling
    isempty(C) || begin
        cl, cu = coupling_bounds(model)
        @constraint(optimization_model, c_lbs, cl .<= C * x) # coupling lower bounds
        @constraint(optimization_model, c_ubs, C * x .<= cu) # coupling upper bounds
    end

    return optimization_model
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
    lower_bound::Maybe{Real} = nothing,
    upper_bound::Maybe{Real} = nothing,
)
    isnothing(lower_bound) || set_normalized_rhs(opt_model[:lbs][vidx], -lower_bound)
    isnothing(upper_bound) || set_normalized_rhs(opt_model[:ubs][vidx], upper_bound)
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

"""
$(TYPEDSIGNATURES)

From the optimized model, returns a vector of values for the selected
`semantics`. If the model did not solve, returns `nothing`.

# Example
```
values_vec(Val(:reaction), model, flux_balance_analysis(model, ...)) # in order of reactions(model)
```
"""
function values_vec(semantics::Val{Semantics}, res::AbstractResult) where {Semantics}
    sem = Accessors.Internal.get_semantics(semantics)
    isnothing(sem) && throw(DomainError(semantics, "Unknown semantics"))
    (_, _, _, sem_varmtx) = sem
    is_solved(res.opt_model) ? sem_varmtx(res.model)' * value.(res.opt_model[:x]) : nothing
end

"""
$(TYPEDSIGNATURES)

Convenience variant of [`values_vec`](@ref).

# Example
```
values_vec(:reaction, model, flux_balance_analysis(model, ...)) # in order of reactions(model)
```
"""
values_vec(semantics::Symbol, res::AbstractResult) = values_vec(Val(semantics), res)

"""
$(TYPEDSIGNATURES)

A pipeable variant of the convenience variant of [`values_vec`](@ref).

# Example
```
flux_balance_analysis(model, ...) |> values_vec(:reaction)
```
"""
values_vec(semantics::Symbol) = res -> values_vec(Val(semantics), res)

"""
$(TYPEDSIGNATURES)

From the optimized model, returns a dictionary mapping semantic IDs to their
solved values for the selected `semantics`. If the model did not solve, returns
`nothing`.

# Example
```
values_dict(Val(:reaction), flux_balance_analysis(model, ...))
```
"""
function values_dict(semantics::Val{Semantics}, res::AbstractResult) where {Semantics}
    sem = Accessors.Internal.get_semantics(semantics)
    isnothing(sem) && throw(DomainError(semantics, "Unknown semantics"))
    (ids, _, _, sem_varmtx) = sem
    is_solved(res.opt_model) ?
    Dict(ids(res.model) .=> sem_varmtx(res.model)' * value.(res.opt_model[:x])) : nothing
end

"""
$(TYPEDSIGNATURES)

Convenience variant of [`values_dict`](@ref).

# Example
```
values_dict(:reaction, model, flux_balance_analysis(model, ...))
```
"""
values_dict(semantics::Symbol, res::AbstractResult) = values_dict(Val(semantics), res)

"""
$(TYPEDSIGNATURES)

A pipeable variant of the convenience variant of [`values_dict`](@ref).

# Example
```
flux_balance_analysis(model, ...) |> values_dict(:reaction)
```
"""
values_dict(semantics::Symbol) = res -> values_dict(Val(semantics), res)

@export_locals
end
