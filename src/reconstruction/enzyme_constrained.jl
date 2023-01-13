"""
$(TYPEDSIGNATURES)

Wrap a model into a [`EnzymeConstrainedModel`](@ref), following the structure
given by GECKO algorithm (see [`EnzymeConstrainedModel`](@ref) documentation for
details). Multiple capacity constraints can be placed on the model using the
kwargs.

# Arguments
- a `model` that implements the accessors `gene_product_molar_mass`,
  `reaction_isozymes`, `gene_product_lower_bound`, `gene_product_upper_bound`.
- `gene_product_mass_group` is a dict that returns a vector of gene IDs
  associated with each a named capacity constraint. By default, all gene
  products belong to group `"uncategorized"`, which is the behavior of original
  GECKO.
- `gene_product_mass_group_bound` is a dict that returns the capacity
  limitation for a given named capacity constraint.

# Example
```
ecmodel = make_enzyme_constrained_model(
    model;
    gene_product_mass_group = Dict(
        "membrane" => ["e1", "e2"],
        "total" => ["e1", "e2", "e3"],
    ),
    gene_product_mass_group_bound = Dict(
        "membrane" => 0.2,
        "total" => 0.5,
    ),
)
```

# Notes
Reactions with no turnover number data, or non-enzymatic reactions that should
be ignored, must have `nothing` in the `gene_associations` field of the
associated reaction.
"""
function make_enzyme_constrained_model(
    model::AbstractMetabolicModel;
    gene_product_mass_group::Dict{String,Vector{String}} = Dict(
        "uncategorized" => genes(model),
    ),
    gene_product_mass_group_bound::Dict{String,Float64} = Dict("uncategorized" => 0.5),
)
    gpb_(gid) = (gene_product_lower_bound(model, gid), gene_product_upper_bound(model, gid))

    gpmm_(gid) = gene_product_molar_mass(model, gid)

    columns = Vector{Types._EnzymeConstrainedReactionColumn}()
    coupling_row_reaction = Int[]
    coupling_row_gene_product = Int[]

    gids = genes(model)
    (lbs, ubs) = bounds(model)
    rids = variables(model)

    gene_name_lookup = Dict(gids .=> 1:length(gids))
    gene_row_lookup = Dict{Int,Int}()

    for i = 1:n_variables(model)
        isozymes = reaction_isozymes(model, rids[i])

        if isnothing(isozymes)
            push!(
                columns,
                Types._EnzymeConstrainedReactionColumn(i, 0, 0, 0, lbs[i], ubs[i], []),
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
                        Types._EnzymeConstrainedReactionColumn(
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
        for (mg, mg_gids) in gene_product_mass_group #  each gid can belong to multiple mass groups
            gid ∉ mg_gids && continue
            if haskey(mg_gid_lookup, mg)
                push!(mg_gid_lookup[mg], gid)
            else
                mg_gid_lookup[mg] = [gid]
            end
        end
    end
    coupling_row_mass_group = Vector{Types._EnzymeConstrainedCapacity}()
    for (grp, gs) in mg_gid_lookup
        idxs = [gene_row_lookup[x] for x in Int.(indexin(gs, gids))]
        mms = gpmm_.(gs)
        push!(
            coupling_row_mass_group,
            Types._EnzymeConstrainedCapacity(
                grp,
                idxs,
                mms,
                gene_product_mass_group_bound[grp],
            ),
        )
    end

    EnzymeConstrainedModel(
        [
            Types.enzyme_constrained_column_reactions(columns, model)' * objective(model)
            spzeros(length(coupling_row_gene_product))
        ],
        columns,
        coupling_row_reaction,
        collect(zip(coupling_row_gene_product, gpb_.(gids[coupling_row_gene_product]))),
        coupling_row_mass_group,
        model,
    )
end
