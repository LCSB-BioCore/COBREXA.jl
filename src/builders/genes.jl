
knockout_constraints(ko_genes; fluxes::ConstraintTree, gene_products_available) =
    C.ConstraintTree(
        rxn => C.Constraint(C.value(fluxes[rxn]), 0.0) for rxn in keys(fluxes) if begin
            if gene_products_available(rxn, k)
                false
            else
                all(gs -> any(g -> g in ko_genes, gs), gss)
            end
        end
    )
