"""
$(TYPEDSIGNATURES)

A helper function that returns the reaction ids that are active. Active reaction
have thermodynamic data AND a flux bigger than `small_flux_tol` AND are not
ignored.
"""
active_reaction_ids(model::MaxMinDrivingForceModel) = filter(
    rid ->
        haskey(model.reaction_standard_gibbs_free_energies, rid) &&
            abs(get(model.flux_solution, rid, model.small_flux_tol / 2)) >
            model.small_flux_tol &&
            !(rid in model.ignore_reaction_ids),
    reaction_ids(model),
)

"""
$(TYPEDSIGNATURES)

Helper function that returns the unmangled variable IDs.
"""
original_variables(model::MaxMinDrivingForceModel) =
    ["mmdf"; metabolite_ids(model); reaction_ids(model)]
