"""
$(TYPEDSIGNATURES)

Run flux balance analysis (FBA) on the `model`, optionally specifying
`modifications` to the problem.  Basically, FBA solves this optimization
problem:
```
max cᵀx
s.t. S x = b
     xₗ ≤ x ≤ xᵤ
```
See "Orth, J., Thiele, I. & Palsson, B. What is flux balance analysis?. Nat
Biotechnol 28, 245-248 (2010). https://doi.org/10.1038/nbt.1614" for more
information.

The `optimizer` must be set to a `JuMP`-compatible optimizer, such as
`GLPK.Optimizer` or `Tulip.Optimizer`.

Optionally, you may specify one or more modifications to be applied to the model
before the analysis, such as [`set_objective_sense`](@ref),
[`set_optimizer`](@ref), [`set_optimizer_attribute`](@ref), and
[`silence`](@ref).

Returns a tree with the optimization solution of the same shape as the model
defined by [`fbc_model_constraints`](@ref).

# Example
```
model = load_model("e_coli_core.json")
solution = flux_balance(model, GLPK.optimizer)
```
"""
function flux_balance(model::A.AbstractFBCModel, optimizer; kwargs...)
    constraints = fbc_model_constraints(model)
    optimized_constraints(
        constraints;
        objective = constraints.objective.value,
        optimizer,
        kwargs...,
    )
end

"""
$(TYPEDSIGNATURES)

Pipe-able overload of [`flux_balance`](@ref).
"""
flux_balance(optimizer; modifications = []) = m -> flux_balance(m, optimizer; modifications)

export flux_balance
