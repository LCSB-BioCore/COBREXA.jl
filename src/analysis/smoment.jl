
"""
    make_smoment_model(
        model::MetabolicModel;
        reaction_isozymes::Function,
        gene_product_molar_mass::Function,
        total_enzyme_capacity::Float64,
    )

Construct a model with a structure given by sMOMENT algorithm; returns a
[`SMomentModel`](@ref) (see the documentation for details.

`reaction_isozymes` parameter is a function that returns a single isozyme for
each reaction, or `nothing` if the reaction is not enzymatic. If the reaction
has multiple isozymes, use [`smoment_isozyme_score`](@ref) to select the "best"
one, as recommended by the sMOMENT approach.

`gene_product_molar_mass` parameter is a function that returns a molar mass of
each gene product (relative to `total_enzyme_capacity` and the specified
kcats), as specified by sMOMENT.

`total_enzyme_capacity` is the maximum "enzyme capacity" consumption of the
model.
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
            continue
        end
        # pick a new row for "arm reaction" coupling
        push!(coupling_row_reaction, i)
        coupling_row = length(coupling_row_reaction)

        mw = sum(
            gene_product_molar_mass(gid) * ps for (gid, ps) in isozyme.gene_product_count
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

        if max(lbs[i], ubs[i]) > 0 && isozyme.kcat_forward > _constants.tolerance
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

    return SMomentModel(columns, coupling_row_reaction, total_enzyme_capacity, model)
end
