
"""
    make_optimization_model(
        model::MetabolicModel,
        optimizer;
        sense = MAX_SENSE,
    )

Convert `MetabolicModel`s to a JuMP model, place objectives and the equality
constraint.

Here coupling means inequality constraints coupling multiple variables together.
"""
function make_optimization_model(model::MetabolicModel, optimizer; sense = MAX_SENSE)

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
    isempty(C) || @constraint(optimization_model, c_lbs, cl .<= C * x) # coupling lower bounds
    isempty(C) || @constraint(optimization_model, c_ubs, C * x .<= cu) # coupling upper bounds

    return optimization_model
end

"""
    is_solved(opt_model)

Return `true` if `opt_model` solved successfully (solution is optimal or locally
optimal).  Return `false` if any other termination status is reached.
Termination status is defined in the documentation of `JuMP`.
"""
is_solved(opt_model) = termination_status(opt_model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]

"""
    optimize_objective(opt_model)::Maybe{Float64}

Shortcut for running JuMP `optimize!` on a model and returning the objective
value, if solved.
"""
function optimize_objective(opt_model)::Maybe{Float64}
    optimize!(opt_model)
    solved_objective_value(opt_model)
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

"""
    solved_objective_value(opt_model)::Maybe{Float64}

Returns the current objective value of a model, if solved.

# Example
```
solved_objective_value(flux_balance_analysis(model, ...))
```
"""
solved_objective_value(opt_model)::Maybe{Float64} =
    is_solved(opt_model) ? objective_value(opt_model) : nothing

"""
    flux_vector(opt_model)::Maybe{Vector{Float64}}

Returns a vector of fluxes of the model, if solved.

# Example
```
flux_vector(flux_balance_analysis(model, ...))
```
"""
flux_vector(model::MetabolicModel, opt_model)::Maybe{Vector{Float64}} =
    is_solved(opt_model) ? reaction_flux(model)' * value.(opt_model[:x]) : nothing

"""
    flux_dict(model::MetabolicModel, opt_model)::Maybe{Dict{String, Float64}, Nothing}

Returns the fluxes of the model as a reaction-keyed dictionary, if solved.

# Example
```
flux_dict(model, flux_balance_analysis(model, ...))
```
"""
flux_dict(model::MetabolicModel, opt_model)::Maybe{Dict{String,Float64}} =
    is_solved(opt_model) ?
    Dict(fluxes(model) .=> reaction_flux(model)' * value.(opt_model[:x])) : nothing

"""
    flux_dict(model::MetabolicModel)

A pipeable variant of `flux_dict`.

# Example
```
flux_balance_analysis(model, ...) |> flux_dict(model)
```
"""
flux_dict(model::MetabolicModel) = opt_model -> flux_dict(model, opt_model)
