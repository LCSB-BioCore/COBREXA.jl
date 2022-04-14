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

        reaction_coupling_row =
            length(isozymes) > 1 ? begin
                push!(coupling_row_reaction, i)
                length(coupling_row_reaction)
            end : 0

        for (iidx, isozyme) in enumerate(isozymes)
            for (lb, ub, kcat, dir) in [
                (-ubs[i], -lbs[i], isozyme.kcat_reverse, -1),
                (lbs[i], ubs[i], isozyme.kcat_forward, 1),
            ]
                if max(lb, ub) > 0 && kcat > _constants.tolerance
                    push!(
                        columns,
                        _gecko_column(
                            i,
                            iidx,
                            dir,
                            reaction_coupling_row,
                            max(lb, 0),
                            ub,
                            _gecko_make_gene_product_coupling(
                                isozyme.gene_product_count,
                                kcat,
                                gene_name_lookup,
                                gene_row_lookup,
                                coupling_row_gene_product,
                            ),
                            _gecko_make_mass_group_coupling(
                                isozyme.gene_product_count,
                                kcat,
                                gene_mass_group,
                                gene_product_mass,
                                mass_group_lookup,
                                coupling_row_mass_group,
                            ),
                        ),
                    )
                end
            end
        end
    end

    coupling_row_mass_group = return GeckoModel(
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

function _gecko_make_mass_group_coupling(
    gene_product_count::Dict{String,Int},
    kcat::Float64,
    gene_row::Function,
    gene_product_mass::Function,
    row_lookup::Dict{String,Int},
    rows::Vector{String},
)
    gp_groups = gene_row.(keys(gene_product_count))
    gp_mass = gene_product_mass.(keys(gene_product_count))
    groups = unique(filter(!isnothing, gp_groups))
    group_idx = Dict(groups .=> 1:length(groups))
    vals = [0.0 for _ in groups]

    for (gpg, mass) in zip(gp_groups, gp_mass)
        if !isnothing(gpg)
            vals[group_idx[gpg]] += mass / kcat
        end
    end

    collect(
        isnothing(group) ? 0 : begin
            if !haskey(row_lookup, group)
                push!(rows, group)
                row_lookup[group] = length(rows)
            end
            (row_lookup[group], val)
        end for (group, val) in zip(groups, vals)
    )
end
