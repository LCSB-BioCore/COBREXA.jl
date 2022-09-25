unwrap_model(model::SMomentModel) = model.inner

"""
$(TYPEDSIGNATURES)

Return a stoichiometry of the [`SMomentModel`](@ref). The enzymatic reactions
are split into unidirectional forward and reverse ones.
"""
stoichiometry(model::SMomentModel) =
    stoichiometry(model.inner) * _smoment_column_reactions(model)

"""
$(TYPEDSIGNATURES)

Reconstruct an objective of the [`SMomentModel`](@ref).
"""
objective(model::SMomentModel) = _smoment_column_reactions(model)' * objective(model.inner)

"""
$(TYPEDSIGNATURES)

Returns the internal reactions in a [`SMomentModel`](@ref) (these may be split
to forward- and reverse-only parts; reactions IDs are mangled accordingly with
suffixes).
"""
reactions(model::SMomentModel) =
    let inner_reactions = reactions(model.inner)
        [
            _smoment_reaction_name(inner_reactions[col.reaction_idx], col.direction) for
            col in model.columns
        ]
    end

"""
$(TYPEDSIGNATURES)

The number of reactions (including split ones) in [`SMomentModel`](@ref).
"""
n_reactions(model::SMomentModel) = length(model.columns)

"""
$(TYPEDSIGNATURES)

Return the variable bounds for [`SMomentModel`](@ref).
"""
bounds(model::SMomentModel) =
    ([col.lb for col in model.columns], [col.ub for col in model.columns])

"""
$(TYPEDSIGNATURES)

Get the mapping of the reaction rates in [`SMomentModel`](@ref) to the original
fluxes in the wrapped model.
"""
reaction_flux(model::SMomentModel) =
    _smoment_column_reactions(model)' * reaction_flux(model.inner)

"""
$(TYPEDSIGNATURES)

Return the coupling of [`SMomentModel`](@ref). That combines the coupling of
the wrapped model, coupling for split reactions, and the coupling for the total
enzyme capacity.
"""
coupling(model::SMomentModel) = vcat(
    coupling(model.inner) * _smoment_column_reactions(model),
    [col.capacity_required for col in model.columns]',
)

"""
$(TYPEDSIGNATURES)

Count the coupling constraints in [`SMomentModel`](@ref) (refer to
[`coupling`](@ref) for details).
"""
n_coupling_constraints(model::SMomentModel) = n_coupling_constraints(model.inner) + 1

"""
$(TYPEDSIGNATURES)

The coupling bounds for [`SMomentModel`](@ref) (refer to [`coupling`](@ref) for
details).
"""
coupling_bounds(model::SMomentModel) =
    let (iclb, icub) = coupling_bounds(model.inner)
        (vcat(iclb, [0.0]), vcat(icub, [model.total_enzyme_capacity]))
    end
