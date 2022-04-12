"""
    protein_dict(model::GeckoModel, opt_model)

Return a dictionary mapping protein concentrations to their ids. The argument
`opt_model` is a solved optimization problem, typically returned by
[`flux_balance_analysis`](@ref).
"""
protein_dict(model::GeckoModel, opt_model) =
    is_solved(opt_model) ?
    Dict(
        model.gene_ids .=> value.(opt_model[:x][(length(model.irrev_reaction_ids)+1):end]),
    ) : nothing

"""
    protein_dict(model::GeckoModel)

A pipe-able variant of `protein_dict`.
"""
protein_dict(model::GeckoModel) = x -> protein_dict(model, x)

"""
    get_genes_with_kcats(rid_isozymes::Dict{String, Vector{Isozyme}})

Return all protein (gene ids) that have a kcat from `model` based on `reaction_kcats` field.
Assume that if a reaction has a kcat then each isozyme has a kcat.
"""
function get_genes_with_kcats(rid_isozymes::Dict{String,Vector{Isozyme}})
    gids = String[]
    for isozymes in values(rid_isozymes)
        for isozyme in isozymes
            append!(gids, collect(keys(isozyme.stoichiometry)))
        end
    end
    return unique(gids)
end

"""
    remove_low_expressed_isozymes!(
        model::StandardModel,
        rid_isozymes =  Dict{String, Vector{Isozyme}}()
        gid_measurements = Dict(),
    )

Modify `rid_isozymes` in place by keeping only the highest expressed isozyme.
"""
function remove_low_expressed_isozymes!(
    model::StandardModel,
    rid_isozymes = Dict{String,Vector{Isozyme}}(),
    gid_measurements = Dict(),
)

    for (rid, isozymes) in rid_isozymes
        measured_proteins = Float64[]
        for isozyme in isozymes
            gid_stoich = isozyme.stoichiometry
            push!(
                measured_proteins,
                sum(
                    map(
                        *,
                        collect(values(gid_stoich)),
                        [get(gid_measurements, gid, 0.0) for gid in keys(gid_stoich)],
                        [model.genes[gid].molar_mass for gid in keys(gid_stoich)],
                    ),
                ),
            )
        end
        idx = argmax(measured_proteins)
        rid_isozymes[rid] = [rid_isozymes[rid][idx]]
    end

    return nothing
end
