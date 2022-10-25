"""
$(TYPEDSIGNATURES)

Wrap a model into a [`GeckoModel`](@ref), following the structure given by
GECKO algorithm (see [`GeckoModel`](@ref) documentation for details).

# Arguments

- `reaction_isozymes` is a function that returns a vector of [`Isozyme`](@ref)s
  for each reaction, or empty vector if the reaction is not enzymatic.
- `gene_product_bounds` is a function that returns lower and upper bound for
  concentration for a given gene product (specified by the same string gene ID as in
  `reaction_isozymes`), as `Tuple{Float64,Float64}`.
- `gene_product_molar_mass` is a function that returns a numeric molar mass of
  a given gene product specified by string gene ID.
- `gene_product_mass_group` is a function that returns a string group identifier for a
  given gene product, again specified by string gene ID. By default, all gene
  products belong to group `"uncategorized"` which is the behavior of original
  GECKO.
- `gene_product_mass_group_bound` is a function that returns the maximum mass for a given
  mass group.

Alternatively, all function arguments may be also supplied as dictionaries that
provide the same data lookup.
"""
function make_gecko_model(
    model::AbstractMetabolicModel;
    reaction_isozymes::Union{Function,Dict{String,Vector{Isozyme}}},
    gene_product_bounds::Union{Function,Dict{String,Tuple{Float64,Float64}}},
    gene_product_molar_mass::Union{Function,Dict{String,Float64}},
    gene_product_mass_group::Union{Function,Dict{String,String}} = _ -> "uncategorized",
    gene_product_mass_group_bound::Union{Function,Dict{String,Float64}},
)
    ris_ =
        reaction_isozymes isa Function ? reaction_isozymes :
        (rid -> get(reaction_isozymes, rid, []))
    gpb_ =
        gene_product_bounds isa Function ? gene_product_bounds :
        (gid -> gene_product_bounds[gid])
    gpmm_ =
        gene_product_molar_mass isa Function ? gene_product_molar_mass :
        (gid -> gene_product_molar_mass[gid])
    gmg_ =
        gene_product_mass_group isa Function ? gene_product_mass_group :
        (gid -> gene_product_mass_group[gid])
    gmgb_ =
        gene_product_mass_group_bound isa Function ? gene_product_mass_group_bound :
        (grp -> gene_product_mass_group_bound[grp])
    # ...it would be nicer to have an overload for this, but kwargs can't be used for dispatch

    columns = Vector{Types._GeckoReactionColumn}()
    coupling_row_reaction = Int[]
    coupling_row_gene_product = Int[]

    gids = genes(model)
    (lbs, ubs) = bounds(model)
    rids = reactions(model)

    gene_name_lookup = Dict(gids .=> 1:length(gids))
    gene_row_lookup = Dict{Int,Int}()

    for i = 1:n_reactions(model)
        isozymes = ris_(rids[i])
        if isempty(isozymes)
            push!(columns, Types._GeckoReactionColumn(i, 0, 0, 0, lbs[i], ubs[i], []))
            continue
        end

        # loop over both directions for all isozymes
        for (lb, ub, kcatf, dir) in [
            (-ubs[i], -lbs[i], i -> i.kcat_reverse, -1),
            (lbs[i], ubs[i], i -> i.kcat_forward, 1),
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
                        end for (gene, stoich) in isozyme.gene_product_count if
                        haskey(gene_name_lookup, gene)
                    )

                    # make a new column
                    push!(
                        columns,
                        Types._GeckoReactionColumn(
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
        mg = gmg_(gid)
        if haskey(mg_gid_lookup, mg)
            push!(mg_gid_lookup[mg], gid)
        else
            mg_gid_lookup[mg] = [gid]
        end
    end
    coupling_row_mass_group = Vector{Types._GeckoCapacity}()
    for (grp, gs) in mg_gid_lookup
        idxs = [gene_row_lookup[x] for x in Int.(indexin(gs, gids))]
        mms = gpmm_.(gs)
        push!(coupling_row_mass_group, Types._GeckoCapacity(grp, idxs, mms, gmgb_(grp)))
    end

    GeckoModel(
        [
            Types.gecko_column_reactions(columns, model)' * objective(model)
            spzeros(length(coupling_row_gene_product))
        ],
        columns,
        coupling_row_reaction,
        collect(zip(coupling_row_gene_product, gpb_.(gids[coupling_row_gene_product]))),
        coupling_row_mass_group,
        model,
    )
end
