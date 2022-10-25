
"""
$(TYPEDSIGNATURES)

A simple helper to pick a single wrapped model
"""
function unwrap_model(a::ModelWrapper)
    missing_impl_error(unwrap_model, (a,))
end

#
# IMPORTANT
#
# The list of inherited functions must be synced with the methods available for [`AbstractMetabolicModel`](@ref).
#

@inherit_model_methods_fn ModelWrapper () unwrap_model () reactions metabolites stoichiometry bounds balance objective fluxes n_fluxes reaction_flux coupling n_coupling_constraints coupling_bounds genes n_genes precache!

@inherit_model_methods_fn ModelWrapper (rid::String,) unwrap_model (rid,) reaction_gene_association reaction_subsystem reaction_stoichiometry reaction_annotations reaction_notes

@inherit_model_methods_fn ModelWrapper (mid::String,) unwrap_model (mid,) metabolite_formula metabolite_charge metabolite_compartment metabolite_annotations metabolite_notes

@inherit_model_methods_fn ModelWrapper (gid::String,) unwrap_model (gid,) gene_annotations gene_notes
