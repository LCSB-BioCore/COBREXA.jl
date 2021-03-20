"""
Flux balance analysis, solves the following problem for the input `model`:
```
max cᵀx
s.t. S x = b
     xₗ ≤ x ≤ xᵤ
```
"""
fluxBalanceAnalysis(model::LM, optimizer) where {LM<:AbstractCobraModel} =
    solveLP(model, optimizer; sense = MOI.MAX_SENSE)
