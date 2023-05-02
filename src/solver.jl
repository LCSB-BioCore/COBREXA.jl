
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

    optimization_model = Model(optimizer)

    # make the variables
    n = variable_count(model)
    @variable(optimization_model, x[1:n])

    # bound the variables
    xl, xu = variable_bounds(model)
    @constraint(optimization_model, lbs, xl .<= x) # lower bounds
    @constraint(optimization_model, ubs, x .<= xu) # upper bounds

    # mark the objective
    let obj = objective(model)
        if obj isa AbstractVector
            # linear objective case
            @objective(optimization_model, sense, obj' * x)
        else
            # quadratic objective case
            @objective(optimization_model, sense, x' * obj * [x; 1])
        end
    end

    # go over the semantics and add bounds if there are any
    for (semname, sem) in Accessors.Internal.get_semantics()
        bounds = sem.bounds(model)
        if isnothing(bounds)
            continue
        elseif typeof(bounds) <: AbstractVector{Float64}
            # equality bounds
            constraints =
                @constraint(optimization_model, sem.mapping_matrix(model) * x .== bounds)
            label = Symbol(semname, :_eqs)
            optimization_model[label] = constraints
            set_name.(c, "$label")
        elseif typeof(bounds) <: Tuple{<:AbstractVector{Float64},<:AbstractVector{Float64}}
            # lower/upper interval bounds
            slb, sub = bounds
            smtx = sem.mapping_matrix(model)
            constraints = @constraint(optimization_model, slb .<= smtx * x)
            label = Symbol(semname, :_lbs)
            optimization_model[label] = constraints
            set_name.(c, "$label")
            constraints = @constraint(optimization_model, smtx * x .<= sub)
            label = Symbol(semname, :_ubs)
            optimization_model[label] = constraints
            set_name.(c, "$label")
        else
            # if the bounds returned something weird, complain loudly.
            throw(
                TypeError(
                    :make_optimization_model,
                    "conversion of $(typeof(model)) bounds",
                    Union{
                        Nothing,
                        <:AbstractVector{Float64},
                        Tuple{<:AbstractVector{Float64},<:AbstractVector{Float64}},
                    },
                    typeof(bounds),
                ),
            )
        end
    end

    # add coupling constraints
    C = coupling(model)
    if !isempty(C)
        cl, cu = coupling_bounds(model)
        @constraint(optimization_model, c_lbs, cl .<= C * x) # coupling lower bounds
        @constraint(optimization_model, c_ubs, C * x .<= cu) # coupling upper bounds
    end

    return optimization_model
    #TODO so well, what about having ModelWithResult right from this point? ;D
end

"""
$(TYPEDSIGNATURES)

Return `true` if `opt_model` solved successfully (solution is optimal or locally
optimal).  Return `false` if any other termination status is reached.
Termination status is defined in the documentation of `JuMP`.
"""
is_solved(opt_model::Model) =
    termination_status(opt_model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]

"""
$(TYPEDSIGNATURES)

Variant of is_solved that works with [`ModelWithResult`](@ref).
"""
is_solved(r::ModelWithResult{<:Model}) = is_solved(r.result)

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
"""
solved_objective_value(opt_model::Model)::Maybe{Float64} =
    is_solved(opt_model) ? objective_value(opt_model) : nothing

"""
$(TYPEDSIGNATURES)

Pipeable variant of [`solved_objective_value`](@ref).

# Example
```
flux_balance_analysis(model, ...) |> solved_objective_value
```
"""
solved_objective_value(x::ModelWithResult{<:Model}) = solved_objective_value(x.result)

"""
$(TYPEDSIGNATURES)

Return a vector of all variable values from the solved model, in the same order
given by [`variable_ids`](@ref).

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
    s = Accessors.Internal.semantics(semantics)
    is_solved(res.result) ? s.mapping_matrix(res.model)' * value.(res.result[:x]) : nothing
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
    is_solved(res.result) ? Dict(variable_ids(res.model) .=> value.(res.result[:x])) :
    nothing
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
    s = Accessors.Internal.semantics(semantics)
    is_solved(res.result) ?
    Dict(s.ids(res.model) .=> s.mapping_matrix(res.model)' * value.(res.result[:x])) :
    nothing
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
