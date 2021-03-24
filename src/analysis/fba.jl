"""
    fluxBalanceAnalysis(model::M, optimizer) where {M<:MetabolicModel}

Flux balance analysis solves the following problem for the input `model`:
```
max cᵀx
s.t. S x = b
     xₗ ≤ x ≤ xᵤ
```

Returns a solved model from [`optimizeModel`](@ref).
"""
fluxBalanceAnalysis(model::M, optimizer) where {M<:MetabolicModel} =
    optimizeModel(model, optimizer; sense = MOI.MAX_SENSE)

"""
    fluxBalanceAnalysisVec(args...)::Maybe{Vector{Float64}}

A variant of FBA that returns a vector of fluxes in the same order as reactions
of the model, if the solution is found.

Arguments are passed to [`fluxBalanceAnalysis`](@ref).
"""
function fluxBalanceAnalysisVec(args...)::Maybe{Vector{Float64}}
    (optmodel, vars) = fluxBalanceAnalysis(args...)
    if termination_status(optmodel) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        value.(vars)
    else
        nothing
    end
end

"""
    fluxBalanceAnalysisDict(model::M, args...)::Maybe{Dict{String, Float64}} where {M <: MetabolicModel}

A variant of FBA that returns a dictionary assigning fluxes to reactions, if the solution is found.

Arguments are passed to [`fluxBalanceAnalysis`](@ref).
"""
function fluxBalanceAnalysisDict(
    model::M,
    args...,
)::Maybe{Dict{String,Float64}} where {M<:MetabolicModel}
    v = fluxBalanceAnalysisVec(model, args...)
    if isnothing(v)
        nothing
    else
        Dict(zip(reactions(model), v))
    end
end
