
"""
$(TYPEDSIGNATURES)

Internal helper for systematically naming reactions in [`SimplifiedEnzymeConstrainedModel`](@ref).
"""
simplified_enzyme_constrained_reaction_name(original_name::String, direction::Int) =
    direction == 0 ? original_name :
    direction > 0 ? "$original_name#forward" : "$original_name#reverse"

"""
$(TYPEDSIGNATURES)

Retrieve a utility mapping between reactions and split reactions; rows
correspond to "original" reactions, columns correspond to "split" reactions.
"""
simplified_enzyme_constrained_column_reactions(model::SimplifiedEnzymeConstrainedModel) =
    sparse(
        [col.reaction_idx for col in model.columns],
        1:length(model.columns),
        [col.direction >= 0 ? 1 : -1 for col in model.columns],
        n_variables(model.inner),
        length(model.columns),
    )
