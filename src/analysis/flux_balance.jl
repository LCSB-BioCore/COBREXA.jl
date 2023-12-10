"""
$(TYPEDSIGNATURES)

Make an JuMP model out of `constraints` using [`optimization_model`](@ref)
(most arguments are forwarded there), then apply the `modifications`, optimize
the model, and return either `nothing` if the optimization failed, or `output`
substituted with the solved values (`output` defaults to `constraints`.

For a "nice" version for simpler finding of metabolic model optima, use
[`flux_balance`](@ref).
"""
function optimized_constraints(
    constraints::C.ConstraintTreeElem;
    modifications = [],
    output = constraints,
    kwargs...,
)
    om = optimization_model(constraints; kwargs...)
    for m in modifications
        m(om)
    end
    J.optimize!(om)
    is_solved(om) ? C.constraint_values(output, J.value.(om[:x])) : nothing
end

export optimized_constraints

"""
$(TYPEDSIGNATURES)

Compute an optimal objective-optimizing solution of the given `model`.

Most arguments are forwarded to [`optimized_constraints`](@ref).

Returns a tree with the optimization solution of the same shape as
given by [`fbc_model_constraints`](@ref).
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
