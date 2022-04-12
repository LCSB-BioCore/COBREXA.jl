
"""
    _smoment_reaction_name(original_name::String, direction::Int)

Internal helper for systematically naming reactions in [`SMomentModel`](@ref).
"""
_smoment_reaction_name(original_name::String, direction::Int) =
    direction == 0 ? original_name :
    direction > 0 ? "$original_name#forward" : "$original_name#reverse"

"""
    _smoment_column_reactions(model::SMomentModel)

Retrieve a utility mapping between reactions and split reactions; rows
correspond to "original" reactions, columns correspond to "split" reactions.
"""
_smoment_column_reactions(model::SMomentModel) = sparse(
    [col.reaction_id for col in model.columns],
    1:length(model.columns),
    [col.direction >= 0 ? 1 : -1 for col in model.columns],
    n_reactions(model.inner),
    length(model.columns),
)

"""
    _smoment_reaction_coupling(model::SMomentModel)

Compute the part of the coupling for [`SMomentModel`](@ref) that limits the
"arm" reactions (which group the individual split unidirectional reactions).
"""
_smoment_reaction_coupling(model::SMomentModel) = sparse(
    [col.coupling_row for col in model.columns if col.direction != 0],
    [i for (i, col) in enumerate(model.columns) if col.direction != 0],
    [col.direction for col in model.columns if col.direction != 0],
    _smoment_n_reaction_couplings(model),
    length(model.columns),
)

"""
    _smoment_n_reaction_couplings(model::SMomentModel)

Internal helper for determining the number of required couplings to account for
"arm" reactions.
"""
_smoment_n_reaction_couplings(model::SMomentModel) =
    length(model.coupling_row_reaction)

"""
    _smoment_reaction_coupling_bounds(model::SMomentModel)

Return bounds that limit the "arm" reactions in [`SMomentModel`](@ref). The
values are taken from the "original" inner model.
"""
_smoment_reaction_coupling_bounds(model::SMomentModel) =
    let (lbs, ubs) = bounds(model.inner)
        (lbs[coupling_row_reaction], ubs[coupling_row_reaction])
    end
