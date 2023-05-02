# TODO move this to wrappers
"""
$(TYPEDSIGNATURES)

Wrap a `model` into an [`EnzymeConstrainedModel`](@ref), following the structure
given by the GECKO algorithm (see [`EnzymeConstrainedModel`](@ref) documentation
for details). Multiple mass constraint groups can be placed on the model using
the keyword arguments.

Parameters `gene_product_mass_groups` and `gene_product_mass_group_bounds` specify
the groups of gene products, and the respective total mass limit for each group.
Gene products that are not listed in any gene product mass group are ignored.

For simplicity, specifying the `total_gene_product_mass_bound` argument
overrides the above arguments by internally specifying a single group called
`uncategorized` of all gene products, and acts like the maximum "enzyme
capacity" in the model.

# Example
```
ecmodel = make_enzyme_constrained_model(
    model;
    gene_product_mass_groups = Dict(
        "membrane" => ["e1", "e2"],
        "total" => ["e1", "e2", "e3"],
    ),
    gene_product_mass_group_bounds = Dict(
        "membrane" => 0.2,
        "total" => 0.5,
    ),
)

ecmodel2 = make_enzyme_constrained_model(
    model;
    total_gene_product_mass_bound = 0.5
)
```
"""
function make_enzyme_constrained_model(
    model::AbstractMetabolicModel;
    gene_product_mass_groups::Maybe{Dict{String,Vector{String}}} = nothing,
    gene_product_mass_group_bounds::Maybe{Dict{String,Float64}} = nothing,
    total_gene_product_mass_bound::Maybe{Float64} = nothing,
)
    if !isnothing(total_gene_product_mass_bound)
        gene_product_mass_groups = Dict("uncategorized" => genes(model))
        gene_product_mass_group_bounds =
            Dict("uncategorized" => total_gene_product_mass_bound)
    end
    isnothing(gene_product_mass_groups) &&
        throw(ArgumentError("missing mass group specification"))
    isnothing(gene_product_mass_group_bounds) &&
        throw(ArgumentError("missing mass group bounds"))

    gpb_(gid) = (gene_product_lower_bound(model, gid), gene_product_upper_bound(model, gid))

    gpmm_(gid) = gene_product_molar_mass(model, gid)

    columns = Vector{Wrappers.Internal.EnzymeConstrainedReactionColumn}()
    coupling_row_reaction = Int[]
    coupling_row_gene_product = Int[]

    gids = genes(model)
    (lbs, ubs) = bounds(model)
    rids = variable_ids(model)

    gene_name_lookup = Dict(gids .=> 1:length(gids))
    gene_row_lookup = Dict{Int,Int}()

    for i = 1:variable_count(model)
        isozymes = reaction_isozymes(model, rids[i])

        if isnothing(isozymes)
            push!(
                columns,
                Wrappers.Internal.EnzymeConstrainedReactionColumn(
                    i,
                    0,
                    0,
                    0,
                    lbs[i],
                    ubs[i],
                    [],
                ),
            )
            continue
        end

        # loop over both directions for all isozymes
        for (lb, ub, kcatf, dir) in [
            (-ubs[i], -lbs[i], x -> x.kcat_backward, -1),
            (lbs[i], ubs[i], x -> x.kcat_forward, 1),
        ]
            # The coefficients in the coupling matrix will be categorized in
            # separate rows for negative and positive reactions. Surprisingly,
            # we do not need to explicitly remember the bounds, because the
            # ones taken from the original model are perfectly okay -- the
            # "reverse" direction is unreachable because of individual
            # bounds on split reactions, and the "forward" direction is
            # properly negated in the reverse case to work nicely with the
            # global lower bound.
            reaction_coupling_row =
                ub > 0 && length(isozymes) > 1 ? begin
                    push!(coupling_row_reaction, i)
                    length(coupling_row_reaction)
                end : 0

            # all isozymes in this direction
            for (iidx, isozyme) in enumerate(isozymes)
                kcat = kcatf(isozyme)
                if ub > 0 && kcat > constants.tolerance
                    # prepare the coupling with gene product molar
                    gene_product_coupling = collect(
                        begin
                            gidx = gene_name_lookup[gene]
                            row_idx = if haskey(gene_row_lookup, gidx)
                                gene_row_lookup[gidx]
                            else
                                push!(coupling_row_gene_product, gidx)
                                gene_row_lookup[gidx] =
                                    length(coupling_row_gene_product)
                            end
                            (row_idx, stoich / kcat)
                        end for (gene, stoich) in isozyme.gene_product_stoichiometry if
                        haskey(gene_name_lookup, gene)
                    )

                    # make a new column
                    push!(
                        columns,
                        Wrappers.Internal.EnzymeConstrainedReactionColumn(
                            i,
                            iidx,
                            dir,
                            reaction_coupling_row,
                            max(lb, 0),
                            ub,
                            gene_product_coupling,
                        ),
                    )
                end
            end
        end
    end

    # prepare enzyme capacity constraints
    mg_gid_lookup = Dict{String,Vector{String}}()
    for gid in gids[coupling_row_gene_product]
        for (mg, mg_gids) in gene_product_mass_groups #  each gid can belong to multiple mass groups
            gid ∉ mg_gids && continue
            if haskey(mg_gid_lookup, mg)
                push!(mg_gid_lookup[mg], gid)
            else
                mg_gid_lookup[mg] = [gid]
            end
        end
    end
    coupling_row_mass_group = Vector{Wrappers.Internal.EnzymeConstrainedCapacity}()
    for (grp, gs) in mg_gid_lookup
        idxs = [gene_row_lookup[x] for x in Int.(indexin(gs, gids))]
        mms = gpmm_.(gs)
        push!(
            coupling_row_mass_group,
            Wrappers.Internal.EnzymeConstrainedCapacity(
                grp,
                idxs,
                mms,
                gene_product_mass_group_bounds[grp],
            ),
        )
    end

    EnzymeConstrainedModel(
        [
            Wrappers.Internal.enzyme_constrained_column_reactions(columns, model)' *
            objective(model)
            spzeros(length(coupling_row_gene_product))
        ],
        columns,
        coupling_row_reaction,
        collect(zip(coupling_row_gene_product, gpb_.(gids[coupling_row_gene_product]))),
        coupling_row_mass_group,
        model,
    )
end
