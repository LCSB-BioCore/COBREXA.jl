
"""
    struct _smoment_column

A helper type that describes the contents of [`SMomentModel`](@ref)s.
"""
struct _smoment_column
    reaction_id::Int # number of the corresponding reaction in the inner model
    direction::Int # 0 if "as is" and unique, -1 if reverse-only part, 1 if forward-only part
    coupling_row::Int # number of row in the coupling (0 if direction==0)
    lb::Float64 # must be 0 if the reaction is unidirectional (if direction!=0)
    ub::Float64
    capacity_required::Float64 # must be 0 for bidirectional reactions (if direction==0)
end

# TODO fix the docstring
"""
    mutable struct SMomentModel <: ModelWrapper

Construct an enzyme capacity constrained model see `Bekiaris, Pavlos Stephanos,
and Steffen Klamt. "Automatic construction of metabolic models with enzyme
constraints." BMC bioinformatics, 2020.` for implementation details.

Note, `"§"` is reserved for internal use as a delimiter, no reaction id should
contain that character. Also note, SMOMENT assumes that each reaction only has a
single enzyme (one GRR) associated with it. It is required that a model be
modified to ensure that this condition is met. For ease-of-use,
[`remove_slow_isozymes!`](@ref) is supplied to effect this. Currently only
`modifications` that change attributes of the `optimizer` are supported.
"""
mutable struct SMomentModel <: ModelWrapper
    columns::Vector{_smoment_column}
    coupling_row_reaction::Vector{Int}
    total_enzyme_capacity::Float64

    inner::MetabolicModel
end

unwrap_model(model::SMomentModel) = model.inner

"""
    stoichiometry(model::SMomentModel)

Return a stoichiometry of the [`SMomentModel`](@ref). The enzymatic reactions
are split into unidirectional forward and reverse ones.
"""
stoichiometry(model::SMomentModel) =
    stoichiometry(model.inner) * _smoment_column_reactions(model)

"""
    objective(model::SMomentModel)

Reconstruct an objective of the [`SMomentModel`](@ref).
"""
objective(model::SMomentModel) = _smoment_column_reactions(model)' * objective(model.inner)

"""
    irreversible_reactions(model::SMomentModel)

Returns the irreversible reactions in `model`.
"""
reactions(model::SMomentModel) =
    let inner_reactions = reactions(model.inner)
        [
            _smoment_reaction_name(inner_reactions[col.reaction_id], col.direction) for
            col in model.columns
        ]
    end

"""
    n_reactions(model::SMomentModel)

The number of reactions (including split ones) in [`SMomentModel`](@ref).
"""
n_reactions(model::SMomentModel) = length(model.columns)

"""
    bounds(model::SMomentModel)

Return the variable bounds for [`SMomentModel`](@ref).
"""
bounds(model::SMomentModel) =
    ([col.lb for col in model.columns], [col.ub for col in model.columns])

"""
    reaction_flux(model::SMomentModel)

Get the mapping of the reaction rates in [`SMomentModel`](@ref) to the original
fluxes in the wrapped model.
"""
reaction_flux(model::SMomentModel) =
    reaction_flux(model.inner)' * _smoment_column_reactions(model)

"""
    coupling(model::SMomentModel)

Return the coupling of [`SMomentModel`](@ref). That combines the coupling of
the wrapped model, coupling for split reactions, and the coupling for the total
enzyme capacity.
"""
coupling(model::SMomentModel) = vcat(
    coupling(model.inner),
    _smoment_reaction_coupling(model),
    [col.capacity_required for col in model.columns]',
)

"""
    n_coupling_constraints(model::SMomentModel)

Count the coupling constraints in [`SMomentModel`](@ref) (refer to
[`coupling`](@ref) for details).
"""
n_coupling_constraints(model::SMomentModel) =
    n_coupling_constraints(model.inner) + _smoment_n_reaction_couplings(model) + 1

"""
    coupling_bounds(model::SMomentModel)

The coupling bounds for [`SMomentModel`](@ref) (refer to [`coupling`](@ref) for
details).
"""
coupling_bounds(model::SMomentModel) =
    let
        (ilb, iub) =
            n_coupling_constraints(model.inner), (rlb, rub) =
                _smoment_reaction_coupling_bounds(model)
        (vcat(ilb, rlb, 0), vcat(iub, rub, model.total_enzyme_capacity))
    end
