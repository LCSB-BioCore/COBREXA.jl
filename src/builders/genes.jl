
#TODO maybe separate this into simple boolean eval function and actual builder
knockout_constraint(ko_genes; fluxes::SolutionTree, reaction_gene_association) =
    C.ConstraintTree(
        rxn => C.Constraint(C.value(fluxes[rxn]), 0.0)
        for rxn=keys(fluxes)
        if begin
            gss=reaction_gene_association(rxn)
            if isnothing(gss)
                false
            else
                all(gs -> any(g-> g in ko_genes, gs), gss)
            end
        end
    )
