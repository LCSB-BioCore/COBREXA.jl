
struct _gecko_column
    reaction_idx::Int
    isozyme_idx::Int
    direction::Int
    reaction_coupling_row::Int
    lb::Float64
    ub::Float64
    gene_product_coupling::Vector{Tuple{Int,Float64}}
    mass_group_row::Int
    mass_required::Float64
end

struct GeckoModel <: ModelWrapper
    columns::Vector{_gecko_column}
    coupling_row_reaction::Vector{Int}
    coupling_row_gene_product::Vector{Tuple{Int,Float64}}
    coupling_row_mass_group::Vector{Tuple{String,Float64}} #TODO add to matrices

    inner::MetabolicModel
end

unwrap_model(model::GeckoModel) = model.inner

"""
    stoichiometry(model::GeckoModel)

Return a stoichiometry of the [`GeckoModel`](@ref). The enzymatic reactions are
split into unidirectional forward and reverse ones, each of which may have
multiple variants per isozyme.
"""
stoichiometry(model::GeckoModel) =
    stoichiometry(model.inner) * _gecko_column_reactions(model)

"""
    objective(model::GeckoModel)

Reconstruct an objective of the [`GeckoModel`](@ref), following the objective
of the inner model.
"""
objective(model::GeckoModel) = _gecko_column_reactions(model)' * objective(model.inner)

"""
    reactions(model::GeckoModel)

Returns the internal reactions in a [`GeckoModel`](@ref) (these may be split
to forward- and reverse-only parts with different isozyme indexes; reactions
IDs are mangled accordingly with suffixes).
"""
reactions(model::GeckoModel) =
    let inner_reactions = reactions(model.inner)
        [
            _gecko_reaction_name(
                inner_reactions[col.reaction_idx],
                col.direction,
                col.isozyme_idx,
            ) for col in model.columns
        ]
    end

"""
    reactions(model::GeckoModel)

Returns the number of all irreversible reactions in `model`.
"""
n_reactions(model::GeckoModel) = length(model.columns)

"""
    bounds(model::GeckoModel)

Return variable bounds for [`GeckoModel`](@ref).
"""
bounds(model::GeckoModel) =
    ([col.lb for col in model.columns], [col.ub for col in model.columns])

"""
    reaction_flux(model::GeckoModel)

Get the mapping of the reaction rates in [`GeckoModel`](@ref) to the original
fluxes in the wrapped model.
"""
reaction_flux(model::GeckoModel) =
    _gecko_column_reactions(model)' * reaction_flux(model.inner)

"""
    coupling(model::GeckoModel)

Return the coupling of [`GeckoModel`](@ref). That combines the coupling of
the wrapped model, coupling for split reactions, and the coupling for the total
enzyme capacity.
"""
coupling(model::GeckoModel) = vcat(
    coupling(model.inner) * _gecko_column_reactions(model),
    _gecko_reaction_coupling(model),
    _gecko_gene_product_coupling(model),
)

"""
    n_coupling_constraints(model::GeckoModel)

Count the coupling constraints in [`GeckoModel`](@ref) (refer to
[`coupling`](@ref) for details).
"""
n_coupling_constraints(model::GeckoModel) =
    n_coupling_constraints(model.inner) +
    length(model.coupling_row_reaction) +
    length(model.coupling_row_gene_product)

"""
    coupling_bounds(model::GeckoModel)

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
            [0.0 for _ in model.coupling_row_gene_product],
        ),
        vcat(
            icub,
            rub[model.coupling_row_reaction],
            [c for (i, c) in model.coupling_row_gene_product],
        ),
    )
end

"""
    reaction_flux(model::GeckoModel)

Helper function to get fluxes from optimization problem.
"""
function reaction_flux(model::GeckoModel)
    R = spzeros(n_fluxes(model), n_genes(model) + n_reactions(model))
    for (i, rid) in enumerate(fluxes(model))
        for_idx = findfirst(
            x -> x == rid * "§ARM§FOR" || x == rid * "§FOR",
            model.irrev_reaction_ids,
        )
        rev_idx = findfirst(
            x -> x == rid * "§ARM§REV" || x == rid * "§REV",
            model.irrev_reaction_ids,
        )
        !isnothing(for_idx) && (R[i, for_idx] = 1.0)
        !isnothing(rev_idx) && (R[i, rev_idx] = -1.0)
    end
    return R'
end
