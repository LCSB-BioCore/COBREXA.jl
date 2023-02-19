
"""
    module Solver

Interface of COBREXA to JuMP solvers; mainly recreation of the
`AbstractMetabolicModel`s as JuMP optimization models, and retrieval of solved
values.

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
    #TODO what about ModelWithResult right from this point? ;D
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
solved_objective_value(opt_model::Model)::Maybe{Float64} =
    is_solved(opt_model) ? objective_value(opt_model) : nothing

"""
$(TYPEDSIGNATURES)

Pipeable variant of [`solved_objective_value`](@ref).
"""
solved_objective_value(x::ModelWithResult{<:Model}) =
    ModelWithResult(x.model, solved_objective_value(x.result))

"""
$(TYPEDSIGNATURES)

Return a vector of all variable values from the solved model, in the same order
given by [`variables`](@ref).

# Example
```
flux_balance_analysis(model, ...) |> values_vec
```
"""
function values_vec(res::ModelWithResult{<:Model})
    is_solved(res.result) ? value.(res.result[:x]) : nothing
end

"""
$(TYPEDSIGNATURES)

Return a vector of all semantic variable values in the model, in the order
given by the corresponding semantics.

# Example
```
values_vec(:reaction, flux_balance_analysis(model, ...))
```
"""
function values_vec(semantics::Symbol, res::ModelWithResult{<:Model})
    (_, _, _, sem_varmtx) = Accessors.Internal.semantics(semantics)
    is_solved(res.result) ? sem_varmtx(res.model)' * value.(res.result[:x]) : nothing
end

"""
$(TYPEDSIGNATURES)

A pipeable variant of [`values_vec`](@ref).

# Example
```
flux_balance_analysis(model, ...) |> values_vec(:reaction)
```
"""
values_vec(semantics::Symbol) =
    (res::ModelWithResult{<:Model}) -> values_vec(semantics, res)

"""
$(TYPEDSIGNATURES)

Return a dictionary of all variable values from the solved model mapped
to their IDs.

# Example
```
flux_balance_analysis(model, ...) |> values_dict
```
"""
function values_dict(res::ModelWithResult{<:Model})
    is_solved(res.result) ? Dict(variables(res.model) .=> value.(res.result[:x])) : nothing
end

"""
$(TYPEDSIGNATURES)

From the optimized model, returns a dictionary mapping semantic IDs to their
solved values for the selected `semantics`. If the model did not solve, returns
`nothing`.

# Example
```
values_dict(:reaction, flux_balance_analysis(model, ...))
```
"""
function values_dict(semantics::Symbol, res::ModelWithResult{<:Model})
    (ids, _, _, sem_varmtx) = Accessors.Internal.semantics(semantics)
    is_solved(res.result) ?
    Dict(ids(res.model) .=> sem_varmtx(res.model)' * value.(res.result[:x])) : nothing
end

"""
$(TYPEDSIGNATURES)

A pipeable variant of [`values_dict`](@ref).

# Example
```
flux_balance_analysis(model, ...) |> values_dict(:reaction)
```
"""
values_dict(semantics::Symbol) =
    (res::ModelWithResult{<:Model}) -> values_dict(semantics, res)

@export_locals
end
