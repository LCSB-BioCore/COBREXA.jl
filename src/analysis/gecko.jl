"""
    make_geckomodel(
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
gm = make_geckomodel(
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
function make_geckomodel(
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
        _add_enzyme_variable(
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
