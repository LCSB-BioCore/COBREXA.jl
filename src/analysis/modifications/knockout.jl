"""
    knockout(gene_ids::Vector{String})

A modification that zeroes the bounds of all reactions that would be knocked
out by the specified genes (effectively disables the reactions).
"""
knockout(gene_ids::Vector{String}) =
    (model, optmodel) -> _do_knockout(model, optmodel, gene_ids)

"""
    knockout(gene_id::String)

A helper variant of [`knockout`](@ref) for a single gene.
"""
knockout(gene_id::String) = knockout([gene_id])

"""
    _do_knockout(model::MetabolicModel, opt_model)

Internal helper for knockouts on generic MetabolicModels. This can be
overloaded so that the knockouts may work differently (more efficiently) with
other models.
"""
function _do_knockout(model::MetabolicModel, opt_model, gene_ids::Vector{String})
    for (rxn_num, rxn_id) in enumerate(reactions(model))
        if all([
            any(in.(gene_ids, Ref(conjunction))) for
            conjunction in reaction_gene_association(model, rxn_id)
        ])
            set_bound(rxn_num, opt_model, ub = 0, lb = 0)
        end
    end
end
