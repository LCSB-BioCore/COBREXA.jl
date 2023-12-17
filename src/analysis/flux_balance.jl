
"""
$(TYPEDSIGNATURES)

Compute an optimal objective-optimizing solution of the given `model`.

Most arguments are forwarded to [`optimized_constraints`](@ref).

Returns a tree with the optimization solution of the same shape as
given by [`fbc_model_constraints`](@ref).
"""
function flux_balance_analysis(model::A.AbstractFBCModel, optimizer; kwargs...)
    constraints = build_flux_balance_model(model)
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
flux_balance_analysis(optimizer; modifications = []) = m -> flux_balance_analysis(m, optimizer; modifications)

export flux_balance
