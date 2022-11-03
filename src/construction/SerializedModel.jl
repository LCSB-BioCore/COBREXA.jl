
# this just generates the necessary wrappers

@_serialized_change_unwrap add_reactions
@_serialized_change_unwrap change_bound
@_serialized_change_unwrap change_bounds
@_serialized_change_unwrap remove_metabolite
@_serialized_change_unwrap remove_metabolites
@_serialized_change_unwrap remove_reaction
@_serialized_change_unwrap remove_reactions

"""
$(TYPEDSIGNATURES)

Returns the model stored in the serialized structure.
"""
function unwrap_serialized(model::Serialized)
    precache!(model)
    model.m
end
