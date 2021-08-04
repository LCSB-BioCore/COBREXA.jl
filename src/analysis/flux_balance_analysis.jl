"""
    flux_balance_analysis_vec(args...)::Union{Vector{Float64},Nothing}

A variant of FBA that returns a vector of fluxes in the same order as reactions
of the model, if the solution is found.

Arguments are passed to [`flux_balance_analysis`](@ref).
"""
function flux_balance_analysis_vec(args...; kwargs...)::Union{Vector{Float64},Nothing}
    optmodel = flux_balance_analysis(args...; kwargs...)
    is_solved(optmodel) || return nothing
    value.(optmodel[:x])
end

"""
    flux_balance_analysis_dict(model::MetabolicModel, args...)::Union{Dict{String, Float64},Nothing}

A variant of FBA that returns a dictionary assigning fluxes to reactions, if
the solution is found. Arguments are passed to [`flux_balance_analysis`](@ref).
"""
function flux_balance_analysis_dict(
    model::MetabolicModel,
    args...;
    kwargs...,
)::Union{Dict{String,Float64},Nothing}
    v = flux_balance_analysis_vec(model, args...; kwargs...)
    isnothing(v) && return nothing
    Dict(zip(reactions(model), v))
end

"""
    flux_balance_analysis(
        model::M,
        optimizer;
        modifications = [],
    ) where {M<:MetabolicModel}

Run flux balance analysis (FBA) on the `model` optionally specifying
`modifications` to the problem.  Basically, FBA solves this optimization problem:
```
max cᵀx
s.t. S x = b
     xₗ ≤ x ≤ xᵤ
```
See "Orth, J., Thiele, I. & Palsson, B. What is flux balance analysis?. Nat
Biotechnol 28, 245-248 (2010). https://doi.org/10.1038/nbt.1614" for more
information.

The `optimizer` must be set to a `JuMP`-compatible optimizer, such as
`GLPK.Optimizer` or `Tulip.Optimizer`

Optionally, you may specify one or more modifications to be applied to the
model before the analysis, such as [`change_optimizer_attribute`](@ref),
[`change_objective`](@ref), and [`change_sense`](@ref).

Returns an optimized `JuMP` model.

# Example
```
model = load_model("e_coli_core.json")
solution = flux_balance_analysis(model, GLPK.optimizer)
value.(solution[:x])  # extract flux steady state from the optimizer

biomass_reaction_id = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")

modified_solution = flux_balance_analysis(model, GLPK.optimizer;
    modifications=[change_objective(biomass_reaction_id)])
```

"""
function flux_balance_analysis(
    model::M,
    optimizer;
    modifications = [],
) where {M<:MetabolicModel}

    opt_model = make_optimization_model(model, optimizer)

    for mod in modifications
        mod(model, opt_model)
    end

    optimize!(opt_model)
    return opt_model
end
