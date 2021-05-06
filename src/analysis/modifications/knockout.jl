"""
    knockout(gene_ids::Vector{String})

A modification that zeroes the bounds of all reactions that would be knocked
out by the specified genes (effectively disables the reactions).
"""
function knockout(gene_ids::Vector{String})
    # Dev note: the three nested for loops are inefficient. However:
    # - gene_ids (user input) will be probably only very few items
    # - model.genes[gene_id].reactions are just a few reactions (most genes
    #   don't code for a lot of reactions)
    # - reaction.grr also should only hold few items (reactions aren't coded by
    #   many different combinations of genes)
    # Let's avoid premature optimization for now and see if anyone ever has
    # problems with this.
    return (model, opt_model) -> begin
        if typeof(model) == StandardModel
            all_reactions = reactions(model)
            for gene_id in gene_ids
                for reaction_id in gene_associated_reactions(model, gene_id)
                    if all([
                        any(in.(gene_ids, Ref(gene_array))) for
                        gene_array in reaction_gene_association(model, reaction_id)
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
        else # fallback knockout
            for (rxn_num, rxn_id) in enumerate(reactions(model))
                if all([
                    any(in.(gene_and, Ref(gene_ids))) for
                    gene_and in reaction_gene_association(model, rxn_id)
                ])
                    set_bound(rxn_num, opt_model, ub = 0, lb = 0)
                end
            end
        end
    end
end

"""
    knockout(gene_id::String)

A helper variant of [`knockout`](@ref) for a single gene.
"""
function knockout(gene_id::String)
    return knockout([gene_id])
end
