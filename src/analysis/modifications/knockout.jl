"""
    knockout(gene_ids::Vector{String})

A modification that zeroes the bounds of all reactions that would be knocked
out by the specified genes (effectively disables the reactions).
"""
knockout(gene_ids::Vector{String}) = (model, optmodel) -> _do_knockout(model, optmodel)

"""
    knockout(gene_id::String)

A helper variant of [`knockout`](@ref) for a single gene.
"""
knockout(gene_id::String) = knockout([gene_id])

"""
    _do_knockout(model::MetabolicModel, opt_model)

Internal helper for knockouts on generic MetabolicModels.
"""
function _do_knockout(model::MetabolicModel, opt_model)
    for (rxn_num, rxn_id) in enumerate(reactions(model))
        if all([
            any(in.(gene_ids, Ref(conjunction))) for
            conjunction in reaction_gene_association(model, rxn_id)
        ])
            set_bound(rxn_num, opt_model, ub = 0, lb = 0)
        end
    end
end

"""
    _do_knockout(model::StandardModel, opt_model)

Internal helper for knockouts specialized to StandardModel (uses a prepared index).
"""
function _do_knockout(model::StandardModel, opt_model)
    all_reactions = reactions(model)
    for gene_id in gene_ids
        for reaction_id in model.genes[gene_id].associated_reactions
            if all([
                any(in.(gene_ids, Ref(conjunction))) for
                conjunction in reaction_gene_association(model, reaction_id)
            ])
                set_bound(
                    first(indexin([reaction_id], all_reactions)),
                    opt_model,
                    ub = 0,
                    lb = 0,
                )
            end
        end
    end
end
