
knockout_constraints(;fluxes::C.ConstraintTree, knockout_test) =
    C.ConstraintTree(
        rxn => C.Constraint(C.value(fluxes[rxn]), 0.0) for rxn in keys(fluxes) if knockout_test(rxn)
    )

export knockout_constraints

fbc_gene_knockout_constraints(;fluxes::C.ConstraintTree, genes, fbc_model::A.AbstractFBCModel) = knockout_constraints(;
    fluxes,
    knockout_test =
        rxn -> !A.reaction_gene_products_available(
            rxn,
            g -> not(g in genes)
        ),
)

export fbc_gene_knockout_constraints
