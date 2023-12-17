
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
