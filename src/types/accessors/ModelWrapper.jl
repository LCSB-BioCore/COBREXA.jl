
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
# The list of inherited functions must be synced with the methods available for
# [`AbstractMetabolicModel`](@ref).
#

@inherit_model_methods_fn AbstractModelWrapper () unwrap_model () variables metabolites stoichiometry bounds balance objective reactions n_reactions reaction_variables coupling n_coupling_constraints coupling_bounds genes n_genes precache!

@inherit_model_methods_fn AbstractModelWrapper (rid::String,) unwrap_model (rid,) reaction_gene_associations reaction_subsystem reaction_stoichiometry reaction_annotations reaction_notes reaction_isozymes

eval_reaction_gene_association(w::AbstractModelWrapper, rid::String; kwargs...) =
    eval_reaction_gene_association(unwrap_model(w), rid; kwargs...)

@inherit_model_methods_fn AbstractModelWrapper (mid::String,) unwrap_model (mid,) metabolite_formula metabolite_charge metabolite_compartment metabolite_annotations metabolite_notes

@inherit_model_methods_fn AbstractModelWrapper (gid::String,) unwrap_model (gid,) gene_annotations gene_notes gene_product_molar_mass gene_product_lower_bound gene_product_upper_bound
