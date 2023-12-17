

"""
$(TYPEDSIGNATURES)

Return true if the reaction denoted by `rxn_id` in `model` is a boundary
reaction, otherwise return false. Checks if on boundary by inspecting the number
of metabolites in the reaction stoichiometry. Boundary reactions have only one
metabolite, e.g. an exchange reaction, or a sink/demand reaction.
"""
is_boundary(model::A.AbstractFBCModel, rxn_id::String) = length(keys(A.reaction_stoichiometry(model, rxn_id))) == 1

export is_boundary
