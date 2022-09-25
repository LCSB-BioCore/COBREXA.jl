unwrap_model(model::GeckoModel) = model.inner

"""
$(TYPEDSIGNATURES)

Return a stoichiometry of the [`GeckoModel`](@ref). The enzymatic reactions are
split into unidirectional forward and reverse ones, each of which may have
multiple variants per isozyme.
"""
function stoichiometry(model::GeckoModel)
    irrevS = stoichiometry(model.inner) * COBREXA._gecko_reaction_column_reactions(model)
    enzS = COBREXA._gecko_gene_product_coupling(model)
    [
        irrevS spzeros(size(irrevS, 1), size(enzS, 1))
        -enzS I(size(enzS, 1))
    ]
end

"""
$(TYPEDSIGNATURES)

Return the objective of the [`GeckoModel`](@ref). Note, the objective is with
respect to the internal variables, i.e. [`reactions(model)`](@ref), which are
the unidirectional reactions and the genes involved in enzymatic reactions that
have kinetic data.
"""
objective(model::GeckoModel) = model.objective

"""
$(TYPEDSIGNATURES)

Returns the internal reactions in a [`GeckoModel`](@ref) (these may be split
to forward- and reverse-only parts with different isozyme indexes; reactions
IDs are mangled accordingly with suffixes).
"""
function reactions(model::GeckoModel)
    inner_reactions = reactions(model.inner)
    mangled_reactions = [
        _gecko_reaction_name(
            inner_reactions[col.reaction_idx],
            col.direction,
            col.isozyme_idx,
        ) for col in model.columns
    ]
    [mangled_reactions; genes(model)]
end

"""
$(TYPEDSIGNATURES)

Returns the number of all irreversible reactions in `model` as well as the
number of gene products that take part in enzymatic reactions.
"""
n_reactions(model::GeckoModel) = length(model.columns) + n_genes(model)

"""
$(TYPEDSIGNATURES)

Return variable bounds for [`GeckoModel`](@ref).
"""
function bounds(model::GeckoModel)
    lbs = [
        [col.lb for col in model.columns]
        [lb for (_, (lb, _)) in model.coupling_row_gene_product]
    ]
    ubs = [
        [col.ub for col in model.columns]
        [ub for (_, (_, ub)) in model.coupling_row_gene_product]
    ]
    (lbs, ubs)
end

"""
$(TYPEDSIGNATURES)

Get the mapping of the reaction rates in [`GeckoModel`](@ref) to the original
fluxes in the wrapped model.
"""
function reaction_flux(model::GeckoModel)
    rxnmat = _gecko_reaction_column_reactions(model)' * reaction_flux(model.inner)
    [
        rxnmat
        spzeros(n_genes(model), size(rxnmat, 2))
    ]
end

"""
$(TYPEDSIGNATURES)

Return the coupling of [`GeckoModel`](@ref). That combines the coupling of the
wrapped model, coupling for split (arm) reactions, and the coupling for the total
enzyme capacity.
"""
function coupling(model::GeckoModel)
    innerC = coupling(model.inner) * _gecko_reaction_column_reactions(model)
    rxnC = _gecko_reaction_coupling(model)
    enzcap = _gecko_mass_group_coupling(model)
    [
        innerC spzeros(size(innerC, 1), n_genes(model))
        rxnC spzeros(size(rxnC, 1), n_genes(model))
        spzeros(length(model.coupling_row_mass_group), length(model.columns)) enzcap
    ]
end

"""
$(TYPEDSIGNATURES)

Count the coupling constraints in [`GeckoModel`](@ref) (refer to
[`coupling`](@ref) for details).
"""
n_coupling_constraints(model::GeckoModel) =
    n_coupling_constraints(model.inner) +
    length(model.coupling_row_reaction) +
    length(model.coupling_row_mass_group)

"""
$(TYPEDSIGNATURES)

The coupling bounds for [`GeckoModel`](@ref) (refer to [`coupling`](@ref) for
details).
"""
function coupling_bounds(model::GeckoModel)
    (iclb, icub) = coupling_bounds(model.inner)
    (ilb, iub) = bounds(model.inner)
    return (
        vcat(
            iclb,
            ilb[model.coupling_row_reaction],
            [0.0 for _ in model.coupling_row_mass_group],
        ),
        vcat(
            icub,
            iub[model.coupling_row_reaction],
            [grp.group_upper_bound for grp in model.coupling_row_mass_group],
        ),
    )
end

"""
$(TYPEDSIGNATURES)

Return the balance of the reactions in the inner model, concatenated with a vector of
zeros representing the enzyme balance of a [`GeckoModel`](@ref).
"""
balance(model::GeckoModel) =
    [balance(model.inner); spzeros(length(model.coupling_row_gene_product))]

"""
$(TYPEDSIGNATURES)

Return the number of genes that have enzymatic constraints associated with them.
"""
n_genes(model::GeckoModel) = length(model.coupling_row_gene_product)

"""
$(TYPEDSIGNATURES)

Return the gene ids of genes that have enzymatic constraints associated with them.
"""
genes(model::GeckoModel) =
    genes(model.inner)[[idx for (idx, _) in model.coupling_row_gene_product]]

"""
$(TYPEDSIGNATURES)

Return the ids of all metabolites, both real and pseudo, for a [`GeckoModel`](@ref).
"""
metabolites(model::GeckoModel) = [metabolites(model.inner); genes(model) .* "#gecko"]

"""
$(TYPEDSIGNATURES)

Return the number of metabolites, both real and pseudo, for a [`GeckoModel`](@ref).
"""
n_metabolites(model::GeckoModel) = n_metabolites(model.inner) + n_genes(model)
