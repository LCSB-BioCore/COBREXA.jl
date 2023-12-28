
"""
$(TYPEDSIGNATURES)
"""
knockout_constraints(; fluxes::C.ConstraintTree, knockout_test) = C.ConstraintTree(
    rxn => C.Constraint(C.value(fluxes[rxn]), C.EqualTo(0.0)) for
    rxn in keys(fluxes) if knockout_test(rxn)
)

"""
$(TYPEDSIGNATURES)
"""
gene_knockouts(;
    fluxes::C.ConstraintTree,
    ko_genes::Vector{String},
    model::A.AbstractFBCModel,
) = knockout_constraints(;
    fluxes,
    knockout_test = rxn -> begin
        maybe_avail = A.reaction_gene_products_available(
            model,
            string(rxn),
            g -> !(g in ko_genes), # not available if knocked out
        )
        isnothing(maybe_avail) ? false : !maybe_avail # negate here because of knockout_constraints
    end,
)

#TODO remove the bang from here, there's no side effect
"""
$(TYPEDSIGNATURES)
"""
knockout!(ctmodel::C.ConstraintTree, ko_genes::Vector{String}, model::A.AbstractFBCModel) =
    ctmodel * :gene_knockouts^gene_knockouts(; fluxes = ctmodel.fluxes, ko_genes, model)

"""
$(TYPEDSIGNATURES)

Pipe-able variant.
"""
knockout!(ko_genes::Vector{String}, model::A.AbstractFBCModel) =
    m -> knockout!(m, ko_genes, model)
