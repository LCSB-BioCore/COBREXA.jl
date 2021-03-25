"""
    fluxBalanceAnalysis(model::LM, optimizer) where {LM<:MetabolicModel}

Flux balance analysis solves the following problem for the input `model`:
```
max cᵀx
s.t. S x = b
     xₗ ≤ x ≤ xᵤ
```

Returns a solved model from [`optimizeModel`](@ref).
"""
fluxBalanceAnalysis(model::LM, optimizer) where {LM<:MetabolicModel} =
    optimizeModel(model, optimizer; sense = MOI.MAX_SENSE)
