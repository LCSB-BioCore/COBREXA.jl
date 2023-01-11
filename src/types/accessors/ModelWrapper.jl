"""
$(TYPEDSIGNATURES)

A simple helper to pick a single wrapped model
"""
function unwrap_model(a::AbstractModelWrapper)
    missing_impl_error(unwrap_model, (a,))
end

#
# IMPORTANT
#
# The list of inherited functions must be synced with the methods available for [`AbstractMetabolicModel`](@ref).
#

# optimizer based accessors

@inherit_model_methods_fn AbstractModelWrapper () unwrap_model () variables n_variables constraints n_constraints stoichiometry bounds balance objective coupling n_coupling_constraints coupling_bounds reactions n_reactions metabolites n_metabolites reaction_variables genes n_genes precache!

# model interface accessors

@inherit_model_methods_fn AbstractModelWrapper () unwrap_model () reactions n_reactions metabolites n_metabolites reaction_variables genes n_genes precache!

@inherit_model_methods_fn AbstractModelWrapper (rid::String,) unwrap_model (rid,) reaction_gene_associations reaction_subsystem reaction_stoichiometry reaction_annotations reaction_notes

@inherit_model_methods_fn AbstractModelWrapper (mid::String,) unwrap_model (mid,) metabolite_formula metabolite_charge metabolite_compartment metabolite_annotations metabolite_notes

@inherit_model_methods_fn AbstractModelWrapper (gid::String,) unwrap_model (gid,) gene_annotations gene_notes
