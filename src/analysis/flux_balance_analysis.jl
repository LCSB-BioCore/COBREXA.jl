"""
    @flux_balance_analysis_vec

# Example

```
julia> flux_balance_analysis_vec(model, GLPK.Optimizer)
3-element Vector{Float64}:
 1.0
 1.0
 1.0
```
"""
macro flux_balance_analysis_vec(model, optimizer)
    model = esc(model)
    optimizer = esc(optimizer)
    return :(flux_balance_analysis_vec($model, $optimizer))
end

"""
    @flux_balance_analysis_vec

# Example
```
@flux_balance_analysis_vec model Tulip.Optimizer begin
    change_objective(biomass)
    modify_constraint(glucose, -8.0, -8.0)
end
```
"""
macro flux_balance_analysis_vec(model, optimizer, kws)
    model = esc(model)
    optimizer = esc(optimizer)
    fs = [m for m in MacroTools.striplines(kws).args]
    mods = Expr(:vect)
    for f in fs
        push!(mods.args, esc(f))
    end
    kwarg = Expr(:kw, :modifications, mods)
    return :(flux_balance_analysis_vec($model, $optimizer; $kwarg))
end

"""
    @flux_balance_analysis_dict

# Example
```
dict = @flux_balance_analysis_dict model Tulip.Optimizer
```
"""
macro flux_balance_analysis_dict(model, optimizer)
    model = esc(model)
    optimizer = esc(optimizer)
    return :(flux_balance_analysis_dict($model, $optimizer))
end

"""
    @flux_balance_analysis_dict

# Example
```
@flux_balance_analysis_dict model Tulip.Optimizer begin
    change_objective(biomass)
    modify_constraint(glucose, -8.0, -8.0)
end
```
"""
macro flux_balance_analysis_dict(model, optimizer, kws)
    model = esc(model)
    optimizer = esc(optimizer)
    fs = [m for m in MacroTools.striplines(kws).args]
    mods = Expr(:vect)
    for f in fs
        push!(mods.args, esc(f))
    end
    kwarg = Expr(:kw, :modifications, mods)
    return :(flux_balance_analysis_dict($model, $optimizer; $kwarg))
end

"""
    flux_balance_analysis_vec(args...)::Union{Vector{Float64},Nothing}

A variant of FBA that returns a vector of fluxes in the same order as reactions
of the model, if the solution is found.
Arguments are passed to [`flux_balance_analysis`](@ref).
"""
function flux_balance_analysis_vec(args...; kwargs...)::Union{Vector{Float64},Nothing}
    optmodel = flux_balance_analysis(args...; kwargs...)

    COBREXA.JuMP.termination_status(optmodel) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED] ||
        return nothing
    value.(optmodel[:x])
end

"""
    flux_balance_analysis_dict(model::M, args...)::Union{Dict{String, Float64},Nothing} where {M <: MetabolicModel}

A variant of FBA that returns a dictionary assigning fluxes to reactions, if
the solution is found. Arguments are passed to [`flux_balance_analysis`](@ref).
"""
function flux_balance_analysis_dict(
    model::M,
    args...;
    kwargs...,
)::Union{Dict{String,Float64},Nothing} where {M<:MetabolicModel}
    v = flux_balance_analysis_vec(model, args...; kwargs...)
    isnothing(v) && return nothing
    Dict(zip(reactions(model), v))
end

"""
    flux_balance_analysis(
        model::M,
        optimizer;
        modifications = [(model, opt_model) -> nothing],
    ) where {M<:MetabolicModel}

Run flux balance analysis (FBA) on the `model` optionally specifying
`modifications` to the problem.

Effectively solves the optimization problem:
```
max cᵀx
s.t. S x = b
     xₗ ≤ x ≤ xᵤ
```

Optionally, you may specify one or more "modifications" to be applied to the models.
[`change_solver_attribute`](@ref),[`change_objective`](@ref), or
[`change_sense`](@ref) for examples of modifications.

The `optimizer` must be set to perform the analysis, any JuMP solver will work.

Returns a solved JuMP model from [`optimize_model`](@ref).

# Example
```
optimizer = GLPK.Optimizer
model = load_model(StandardModel, "e_coli_core.json")
biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
solved_model = fba(model, optimizer; modifications=[change_objective(biomass)])
```

"""
function flux_balance_analysis(
    model::M,
    optimizer;
    modifications = [(model, opt_model) -> nothing],
) where {M<:MetabolicModel}

    opt_model = make_optimization_model(model, optimizer)

    # support for multiple modification, fallback to single one
    if typeof(modifications) <: AbstractVector
        for mod in modifications
            mod(model, opt_model)
        end
    else
        modifications(model, opt_model)
    end

    COBREXA.JuMP.optimize!(opt_model)
    return opt_model
end
