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

Returns a [`C.ValueTree`](@ref).

# Example
```
model = load_model("e_coli_core.json")
solution = flux_balance_analysis(model, GLPK.optimizer)
```
"""
function flux_balance_analysis(model::A.AbstractFBCModel, optimizer; modifications = [])
    ctmodel = fbc_model_constraints(model)
    flux_balance_analysis(ctmodel, optimizer; modifications)
end

"""
$(TYPEDSIGNATURES)

A variant of [`flux_balance_analysis`](@ref) that takes in a
[`C.ConstraintTree`](@ref) as the model to optimize. The objective is inferred
from the field `objective` in `ctmodel`. All other arguments are forwarded.
"""
function flux_balance_analysis(ctmodel::C.ConstraintTree, optimizer; modifications = [])
    opt_model = optimization_model(
        ctmodel;
        objective = ctmodel.objective.value,
        optimizer,
    )

    for mod in modifications
        mod(ctmodel, opt_model)
    end

    J.optimize!(opt_model)
    
    is_solved(opt_model) || return nothing

    C.ValueTree(ctmodel, J.value.(opt_model[:x]))
end

"""
$(TYPEDSIGNATURES)

Pipe-able variant of [`flux_balance_analysis`](@ref).
"""
flux_balance_analysis(optimizer; modifications = []) =
    m -> flux_balance_analysis(m, optimizer; modifications)

export flux_balance_analysis