"""
    protein_dict(model::GeckoModel, opt_model)

Return a dictionary mapping protein concentrations to their ids.
"""
protein_dict(model::GeckoModel, opt_model) =
    is_solved(opt_model) ?
    last(_map_irrev_to_rev_ids(model.geckodata.reaction_map, value.(opt_model[:x]); protein_ids=model.geckodata.protein_ids)) : nothing


"""
    remove_slow_isozymes!(
        model::StandardModel;
        reaction_kcats = Dict(),
        protein_stoichiometry = Dict(),
        protein_masses = Dict(),
    )

Remove all but the fastest isozyme from each reaction in `model`.
Use the largest kcat (for, rev) for these calculations. Modifies all 
the arguments in place.
"""
function remove_slow_isozymes!(
    model::StandardModel;
    reaction_kcats = Dict(),
    reaction_protein_stoichiometry = Dict(),
    protein_masses = Dict(),
)
    for rid in reactions(model)
        if has_reaction_grr(model, rid) && haskey(reaction_kcats, rid)
            kcat_effs = Float64[]
            grrs = reaction_gene_association(model, rid)
            for (i, grr) in enumerate(grrs)
                push!(
                    kcat_effs,
                    dot(
                        reaction_protein_stoichiometry[rid][i],
                        [protein_masses[gid] for gid in grr],
                    ) / maximum(reaction_kcats[rid][i]),
                )
            end
            idx = argmin(kcat_effs)

            model.reactions[rid].grr = [grrs[idx]]
            reaction_kcats[rid] = [reaction_kcats[rid][idx]]
            reaction_protein_stoichiometry[rid] = [reaction_protein_stoichiometry[rid][idx]]
        end
    end

    curated_gids = String[]
    for rid in reactions(model)
        if has_reaction_grr(model, rid)
            for grr in reaction_gene_association(model, rid)
                append!(curated_gids, grr)
            end
        end
    end
    rm_gids = setdiff(genes(model), curated_gids)
    delete!(model.genes, rm_gids) # remove genes that were deleted

    return nothing
end