
"""
$(TYPEDSIGNATURES)

Run minimization of metabolic adjustment (MOMA) on `model` with respect to
`reference_solution`, which is a dictionary of fluxes. MOMA finds the shortest
Euclidian distance between `reference_solution` and `model` with `modifications`:
```
min Σᵢ (xᵢ - flux_refᵢ)²
s.t. S x = b
     xₗ ≤ x ≤ xᵤ
```
Because the problem has a quadratic objective, a QP solver is required. See
"Daniel, Vitkup & Church, Analysis of Optimality in Natural and Perturbed
Metabolic Networks, Proceedings of the National Academy of Sciences, 2002" for
more details.

Returns a [`C.ValueTree`](@ref), or `nothing` if the solution could not be
found.

# Example
```
model = load_model("e_coli_core.json")

```
"""
function minimal_metabolic_adjustment(
    constraints::C.ConstraintTree,
    reference_solution::Dict{String,Float64},
    optimizer;
    modifications = [],
)
    _constraints =
        constraints *
        :momaobjective^squared_sum_error_objective(
            constraints.fluxes,
            Dict(Symbol(k) => float(v) for (k, v) in reference_solution),
        )

    opt_model = optimization_model(
        _constraints;
        objective = _constraints.momaobjective.value,
        optimizer,
        sense = J.MIN_SENSE,
    )

    for mod in modifications
        mod(constraints, opt_model)
    end

    J.optimize!(opt_model)

    is_solved(opt_model) || return nothing

    C.ValueTree(_constraints, J.value.(opt_model[:x]))
end

"""
$(TYPEDSIGNATURES)

Variant that takes an [`A.AbstractFBCModel`](@ref) as input. All other arguments are forwarded.
"""
function minimal_metabolic_adjustment(model::A.AbstractFBCModel, args...; kwargs...)
    constraints = fbc_model_constraints(model)
    minimal_metabolic_adjustment(constraints, args...; kwargs...)
end

export minimal_metabolic_adjustment
