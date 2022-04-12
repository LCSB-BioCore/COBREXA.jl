"""
    mutable struct GeckoModel <: MetabolicModel

A model that incorporates enzyme capacity and kinetic constraints via the GECKO
formulation. See `Sánchez, Benjamín J., et al. "Improving the phenotype
predictions of a yeast genome‐scale metabolic model by incorporating enzymatic
constraints." Molecular systems biology, 2017.` for implementation details.

Note, since the model uses irreversible reactions internally, `"§FOR"` (for the
forward direction) and `"§REV"` (for the reverse direction) is appended to each
reaction internally. Hence, `"§"` is reserved for internal use as a delimiter,
no reaction id should contain this character. 

To actually run GECKO, call [`flux_balance_analysis`](@ref) on a `GeckoModel`.

# Fields
```
reaction_ids::Vector{String}
irrev_reaction_ids::Vector{String}
metabolites::Vector{String}
gene_ids::Vector{String}
c::SparseVec
S::SparseMat
b::SparseVec
xl::SparseVec
xu::SparseVec
C::SparseMat
cl::Vector{Float64}
cu::Vector{Float64}
```
"""
mutable struct GeckoModel <: MetabolicModel
    reaction_ids::Vector{String}
    irrev_reaction_ids::Vector{String}
    metabolites::Vector{String}
    gene_ids::Vector{String}

    # gecko  
    c::SparseVec
    S::SparseMat
    b::SparseVec
    xl::SparseVec
    xu::SparseVec

    # enzyme capacity constraints
    C::SparseMat
    cl::Vector{Float64}
    cu::Vector{Float64}
end

"""
    stoichiometry(model::GeckoModel)

Return stoichiometry matrix that includes enzymes as metabolites.
"""
stoichiometry(model::GeckoModel) = model.S

"""
    balance(model::GeckoModel)

Return stoichiometric balance.
"""
balance(model::GeckoModel) = model.b

"""
    objective(model::GeckoModel)

Return objective of `model`.
"""
objective(model::GeckoModel) = model.c

"""
    reactions(model::GeckoModel)

Returns the reversible reactions in `model`. For 
the irreversible reactions, use [`irreversible_reactions`][@ref].
"""
reactions(model::GeckoModel) = model.reaction_ids

"""
    n_reactions(model::GeckoModel)

Returns the number of reactions in the model.
"""
n_reactions(model::GeckoModel) = length(model.reaction_ids)

"""
    irreversible_reactions(model::GeckoModel)

Returns the irreversible reactions in `model`.
"""
irreversible_reactions(model::GeckoModel) = model.irrev_reaction_ids

"""
    genes(model::GeckoModel)

Returns the genes (proteins) in the order as they appear as variables in the
model.
"""
genes(model::GeckoModel) = model.gene_ids

"""
    n_genes(model::GeckoModel)

Returns the number of genes in the model.
"""
n_genes(model::GeckoModel) = length(model.gene_ids)

"""
    metabolites(model::GeckoModel)

Return the metabolites in `model`.
"""
metabolites(model::GeckoModel) = model.metabolites

"""
    n_metabolites(model::GeckoModel) = 

Return the number of metabolites in `model`.
"""
n_metabolites(model::GeckoModel) = length(metabolites(model))

"""
    bounds(model::GeckoModel)

Return variable bounds for `GeckoModel`.
"""
bounds(model::GeckoModel) = (model.xl, model.xu)

"""
    coupling(model::GeckoModel)

Coupling constraint matrix for a `GeckoModel`.
"""
coupling(model::GeckoModel) = model.C

"""
    coupling_bounds(model::GeckoModel)

Coupling bounds for a `GeckoModel`.
"""
coupling_bounds(model::GeckoModel) = (model.cl, model.cu)

"""
    reaction_flux(model::MetabolicModel)

Helper function to get fluxes from optimization problem.
"""
function reaction_flux(model::GeckoModel)
    R = spzeros(n_reactions(model), n_genes(model) + length(model.irrev_reaction_ids))
    for (i, rid) in enumerate(reactions(model))
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

"""
    GeckoModel(
        model::StandardModel;
        rid_isozymes = Dict{String, Vector{Isozyme}}(),
        enzyme_capacities = [(),],
    ) 

Construct a `GeckoModel` based on `model` using the kinetic data encoded by 
`rid_isozymes`. Enzyme capacity constraints can be added through `enzyme_capacities`,
which is a vector of tuples. In the first position of the tuple is a list of gene ids,
and the second position is mass upperbound of the sum of these gene ids. 

The units of the fluxes and protein concentration depend on those used in
`rid_isozymes` for the kcats and the molar masses encoded in  the genes of
`model`. Currently only `modifications` that change attributes of the
`optimizer` are supported.

# Example
```
gm = GeckoModel(
    model;
    rid_isozymes,
    enzyme_capacities = [(get_genes_with_kcats(rid_isozymes), total_protein_mass)],
)

opt_model = flux_balance_analysis(
    gm,
    Tulip.Optimizer
)

rxn_fluxes = flux_dict(gm, opt_model)
prot_concens = protein_dict(gm, opt_model)
```
"""
function GeckoModel(
    model::StandardModel;
    rid_isozymes = Dict{String,Vector{Isozyme}}(),
    enzyme_capacities = [()],
)
    S, lb_fluxes, ub_fluxes, reaction_map, metabolite_map =
        _build_irreversible_stoichiometric_matrix(model, rid_isozymes)

    #: find all gene products that have kcats associated with them
    gene_ids = get_genes_with_kcats(rid_isozymes)

    #: size of resultant model
    num_reactions = size(S, 2)
    num_genes = length(gene_ids)
    num_metabolites = size(S, 1)
    num_vars = num_reactions + num_genes

    #: equality lhs
    E_components = ( #TODO add size hints if possible
        row_idxs = Vector{Int}(),
        col_idxs = Vector{Int}(),
        coeffs = Vector{Float64}(),
    )

    for (rid, col_idx) in reaction_map
        original_rid = string(split(rid, "§")[1])

        # skip these entries
        contains(rid, "§ARM") && continue
        !haskey(rid_isozymes, original_rid) && continue

        # these entries have kcats
        if contains(rid, "§ISO")
            iso_num = parse(
                Int,
                replace(
                    first(filter(startswith("ISO"), split(rid, "§")[2:end])),
                    "ISO" => "",
                ),
            )
        else # only one enzyme
            iso_num = 1
        end

        # add all entries to column of matrix
        COBREXA._add_enzyme_variable(
            rid_isozymes,
            iso_num, # only one enzyme
            rid,
            original_rid,
            E_components,
            col_idx,
            gene_ids,
        )
    end

    Se = sparse(
        E_components.row_idxs,
        E_components.col_idxs,
        E_components.coeffs,
        num_genes,
        num_reactions,
    )

    stoich_mat = sparse([
        S zeros(num_metabolites, num_genes)
        Se I(num_genes)
    ])

    #: equality rhs
    b = spzeros(num_metabolites + num_genes)

    #: find objective (assume objective is forward)
    obj_idx_orig = first(findnz(objective(model))[1])
    obj_id_orig = reactions(model)[obj_idx_orig]
    obj_id = obj_id_orig * "§FOR"
    c = spzeros(num_vars)
    obj_idx = reaction_map[obj_id]
    c[obj_idx] = 1.0

    #: inequality constraints
    xl = sparse([lb_fluxes; fill(0.0, num_genes)])
    xu = sparse([ub_fluxes; fill(1000.0, num_genes)])

    #: enzyme capacity constraints
    mw_proteins = [model.genes[pid].molar_mass for pid in gene_ids]
    C = spzeros(length(enzyme_capacities), num_vars)
    cl = spzeros(length(enzyme_capacities))
    cu = spzeros(length(enzyme_capacities))

    for (i, enz_cap) in enumerate(enzyme_capacities)
        enz_idxs = indexin(first(enz_cap), gene_ids)
        C[i, num_reactions.+enz_idxs] .= mw_proteins[enz_idxs]
        cu[i] = last(enz_cap)
    end

    return GeckoModel(
        reactions(model),
        _order_id_to_idx_dict(reaction_map),
        _order_id_to_idx_dict(metabolite_map),
        gene_ids,
        c,
        stoich_mat,
        b,
        xl,
        xu,
        C,
        cl,
        cu,
    )
end

"""
    _add_enzyme_variable(
        rid_isozymes,
        iso_num,
        rid,
        original_rid,
        E_components,
        col_idx,
        gene_ids,
    )

Helper function to add an column into the enzyme stoichiometric matrix.
"""
function _add_enzyme_variable(
    rid_isozymes,
    iso_num,
    rid,
    original_rid,
    E_components,
    col_idx,
    gene_ids,
)
    pstoich = rid_isozymes[original_rid][iso_num].stoichiometry
    kcat =
        contains(rid, "§FOR") ? rid_isozymes[original_rid][iso_num].kcats[1] :
        rid_isozymes[original_rid][iso_num].kcats[2]
    for (pid, pst) in pstoich
        push!(E_components.row_idxs, first(indexin([pid], gene_ids)))
        push!(E_components.col_idxs, col_idx)
        push!(E_components.coeffs, -pst / kcat)
    end
end

"""
    change_bound!(model::GeckoModel, id; lower=nothing, upper=nothing)

Change the bound of variable in `model`. Does not change the bound if respective
bound is `nothing`. Note, for `GeckoModel`s, if the model used to construct the
`GeckoModel` has irreversible reactions, then these reactions will be
permanently irreversible in the model, i.e. changing their bounds to make them
reversible will have no effect.
"""
function change_bound!(model::GeckoModel, id; lower = nothing, upper = nothing)
    gene_idx = first(indexin([id], model.gene_ids))

    if isnothing(gene_idx)
        flux_for_idx = findfirst(
            x -> x == id * "§ARM§FOR" || x == id * "§FOR",
            model.irrev_reaction_ids,
        )
        if !isnothing(flux_for_idx)
            if !isnothing(lower)
                if lower <= 0
                    model.xl[flux_for_idx] = 0
                else
                    model.xl[flux_for_idx] = lower
                end
            end
            if !isnothing(upper)
                if upper <= 0
                    model.xu[flux_for_idx] = 0
                else
                    model.xu[flux_for_idx] = upper
                end
            end
        end

        flux_rev_idx = findfirst(
            x -> x == id * "§ARM§REV" || x == id * "§REV",
            model.irrev_reaction_ids,
        )
        if !isnothing(flux_rev_idx)
            if !isnothing(lower)
                if lower >= 0
                    model.xu[flux_rev_idx] = 0
                else
                    model.xu[flux_rev_idx] = -lower
                end
                if !isnothing(upper)
                    if upper >= 0
                        model.xl[flux_rev_idx] = 0
                    else
                        model.xl[flux_rev_idx] = -upper
                    end
                end
            end
        end
    else
        n = length(model.irrev_reaction_ids)
        !isnothing(lower) && (model.xl[n+gene_idx] = lower)
        !isnothing(upper) && (model.xu[n+gene_idx] = upper)
    end

    return nothing
end

"""
    change_bounds!(model::GeckoModel, ids; lower=fill(nothing, length(ids)), upper=fill(nothing, length(ids)))

Change the bounds of multiple variables in `model` simultaneously. See 
[`change_bound`](@ref) for details.
"""
function change_bounds!(
    model::GeckoModel,
    ids;
    lower = fill(nothing, length(ids)),
    upper = fill(nothing, length(ids)),
)
    for (id, lower, upper) in zip(ids, lower, upper)
        change_bound!(model, id; lower = lower, upper = upper)
    end
end
