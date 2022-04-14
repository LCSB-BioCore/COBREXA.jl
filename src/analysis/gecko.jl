function make_gecko_model(
    model::StandardModel;
    reaction_isozymes::Function,
    gene_product_mass::Function,
    gene_product_limit::Function,
    gene_mass_group::Function = _ -> "uncategorized",
    mass_fraction_limit::Function,
)

    columns = Vector{_gecko_column}()
    coupling_row_reaction = Int[]
    coupling_row_gene_product = Int[]
    coupling_row_mass_group = String[]

    gene_name_lookup = Dict(genes(model) .=> 1:n_genes(model))
    gene_row_lookup = Dict{Int,Int}()
    mass_group_lookup = Dict{String,Int}()

    (lbs, ubs) = bounds(model)
    rids = reactions(model)

    for i = 1:n_reactions(model)
        isozymes = reaction_isozymes(rids[i])
        if isempty(isozymes)
            push!(columns, _gecko_column(i, 0, 0, 0, lbs[i], ubs[i], [], []))
            continue
        end

        # if the reaction has multiple isozymes, it needs extra coupling to
        # ensure that the total rate of the reaction doesn't exceed the
        # "global" limit
        reaction_coupling_row =
            length(isozymes) > 1 ? begin
                push!(coupling_row_reaction, i)
                length(coupling_row_reaction)
            end : 0

        for (iidx, isozyme) in enumerate(isozymes)
            # loop over both directions for all isozymes
            for (lb, ub, kcat, dir) in [
                (-ubs[i], -lbs[i], isozyme.kcat_reverse, -1),
                (lbs[i], ubs[i], isozyme.kcat_forward, 1),
            ]
                if max(lb, ub) > 0 && kcat > _constants.tolerance
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
                            (row_idx, 1 / kcat)
                        end for (gene, count) in isozyme.gene_product_count if
                        haskey(gene_name_lookup, gene)
                    )

                    # prepare the coupling with the mass groups
                    gp_groups = gene_mass_group.(keys(isozyme.gene_product_count))
                    gp_mass = gene_product_mass.(keys(isozyme.gene_product_count))
                    groups = unique(filter(!isnothing, gp_groups))
                    group_idx = Dict(groups .=> 1:length(groups))
                    vals = [0.0 for _ in groups]

                    for (gpg, mass) in zip(gp_groups, gp_mass)
                        if !isnothing(gpg)
                            vals[group_idx[gpg]] += mass / kcat
                        end
                    end

                    mass_group_coupling = collect(
                        isnothing(group) ? 0 :
                        begin
                            if !haskey(mass_group_lookup, group)
                                push!(coupling_row_mass_group, group)
                                mass_group_lookup[group] =
                                    length(coupling_row_mass_group)
                            end
                            (mass_group_lookup[group], val)
                        end for (group, val) in zip(groups, vals)
                    )

                    # make a new column
                    push!(
                        columns,
                        _gecko_column(
                            i,
                            iidx,
                            dir,
                            reaction_coupling_row,
                            max(lb, 0),
                            ub,
                            gene_product_coupling,
                            mass_group_coupling,
                        ),
                    )
                end
            end
        end
    end

    GeckoModel(
        columns,
        coupling_row_reaction,
        collect(
            zip(coupling_row_gene_product, gene_product_limit.(coupling_row_gene_product)),
        ),
        collect(
            zip(coupling_row_mass_group, mass_fraction_limit.(coupling_row_mass_group)),
        ),
        model,
    )
end
