function make_gecko_model(
    model::StandardModel;
    reaction_isozymes::Function,
    reaction_isozyme_masses::Function,
    gene_product_limit::Function,
    reaction_mass_group::Function = _ -> "uncategorized",
    mass_faction_limit::Function,
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
            push!(columns, _gecko_column(i, 0, 0, 0, lbs[i], ubs[i], 0, 0, 0, 0))
            continue
        end

        group = reaction_mass_group(rids[i])
        mass_group_row =
            isnothing(group) ? 0 :
            haskey(massgroup_lookup, group) ? mass_group_lookup[group] :
            begin
                push!(coupling_row_mass_group, group)
                mass_group_lookup[group] = length(coupling_row_mass_group)
            end

        push!(coupling_row_reaction, i)
        reaction_coupling_row = length(coupling_row_reaction)

        masses = group > 0 ? reaction_isozyme_masses(rids[i]) : zeros(length(isozymes))

        for (iidx, isozyme) in enumerate(isozymes)
            if min(lbs[i], ubs[i]) < 0 && isozyme.kcat_reverse > _constants.tolerance
                # reaction can run in reverse
                push!(
                    columns,
                    _gecko_column(
                        i,
                        iidx,
                        -1,
                        reaction_coupling_row,
                        max(-ubs[i], 0),
                        -lbs[i],
                        _gecko_make_gene_product_coupling(
                            isozyme.gene_product_count,
                            isozyme.kcat_reverse,
                            gene_name_lookup,
                            gene_row_lookup,
                            coupling_row_gene_product,
                        ),
                        mass_group_row,
                        masses[iidx] / isozyme.kcat_reverse,
                    ),
                )
            end
            if max(lbs[i], ubs[i]) > 0 && isozyme.kcat_forward > _constants.tolerance
                # reaction can run forward
                push!(
                    columns,
                    _gecko_column(
                        i,
                        iidx,
                        1,
                        reaction_coupling_row,
                        max(lbs[i], 0),
                        ubs[i],
                        _gecko_make_gene_product_coupling(
                            isozyme.gene_product_count,
                            isozyme.kcat_forward,
                            gene_name_lookup,
                            gene_row_lookup,
                            coupling_row_gene_product,
                        ),
                        mass_group_row,
                        masses[iidx] / isozyme.kcat_forward,
                    ),
                )
            end
        end
    end

    coupling_row_mass_group =
        collect(zip(coupling_row_mass_group, mass_fraction_limit.(coupling_row_mass_group)))

    coupling_row_gene_product = collect(
        zip(coupling_row_gene_product, gene_product_limit.(coupling_row_gene_product)),
    )

    return GeckoModel(
        columns,
        coupling_row_reaction,
        coupling_row_gene_product,
        coupling_row_mass_group,
        model,
    )
end

_gecko_make_gene_product_coupling(
    gene_product_count::Dict{String,Int},
    kcat::Float64,
    name_lookup::Dict{String,Int},
    row_lookup::Dict{Int,Int},
    rows::Vector{Int},
) = collect(
    begin
        gidx = name_lookup[gene]
        row_idx = if haskey(row_lookup, gidx)
            row_lookup[gidx]
        else
            push!(rows, gidx)
            row_lookup[gidx] = length(rows)
        end
        (row_idx, 1 / kcat)
    end for (gene, count) in gene_product_count if haskey(name_lookup, gene)
)
