
"""
    make_smoment_model(
        model::MetabolicModel;
        reaction_isozymes::Function,
        gene_product_capacity_required::Function,
        total_enzyme_capacity,
    )

Construct an [`SMomentModel`](@ref) model using the inner `model` and a map of
isozymes.
"""
function make_smoment_model(
    model::MetabolicModel;
    reaction_isozymes::Function,
    gene_product_molar_mass::Function,
    total_enzyme_capacity::Float64,
)
    columns = Vector{_smoment_column}()
    coupling_row_reaction = Int[]

    (lbs, ubs) = bounds(model)
    rids = reactions(model)

    for i = 1:n_reactions(model)
        isozyme = reaction_isozymes(rids[i])
        if isnothing(isozyme)
            # non-enzymatic reaction (or a totally ignored one)
            push!(columns, _smoment_column(i, 0, 0, lbs[i], ubs[i], 0))
        else
            # pick a new row for "arm reaction" coupling
            coupling_row = length(coupling_row_reaction) + 1
            push!(coupling_row_reaction, i)

            mw = sum(
                gene_product_molar_mass(gid) * ps for (gid, ps) in isozyme.stoichiometry
            )

            if min(lbs[i], ubs[i]) < 0 && isozyme.kcat_reverse > _constants.tolerance
                # reaction can run in reverse
                push!(
                    columns,
                    _smoment_column(
                        i,
                        -1,
                        coupling_row,
                        max(-ubs[i], 0),
                        -lbs[i],
                        mw / isozyme.kcat_reverse,
                    ),
                )
            end

            if min(lbs[i], ubs[i]) > 0 && isozyme.kcat_forward > _constants.tolerance
                # reaction can run forward
                push!(
                    columns,
                    _smoment_column(
                        i,
                        1,
                        coupling_row,
                        max(lbs[i], 0),
                        ubs[i],
                        mw / isozyme.kcat_forward,
                    ),
                )
            end
        end
    end

    return SMomentModel(columns, coupling_row_reaction, total_enzyme_capacity, model)
end
