"""
$(TYPEDSIGNATURES)

Wrap a `model` with a structure given by sMOMENT algorithm; returns a
[`SimplifiedEnzymeConstrainedModel`](@ref) (see the documentation for details).
The sMOMENT algorithm only uses one [`Isozyme`](@ref) per reaction. If multiple
isozymes are present in `model`, the "fastest" isozyme will be used. This is
determined based on maximum kcat (forward or backward) divided by mass of the
isozyme. The `total_gene_product_mass_bound` is the maximum "enzyme capacity" in
the model.

# Example
```
ecmodel = make_simplified_enzyme_constrained_model(
    model;
    total_gene_product_mass_bound = 0.5
)
"""
function make_simplified_enzyme_constrained_model(
    model::AbstractMetabolicModel;
    total_gene_product_mass_bound::Float64 = 0.5,
)
    # helper function to rank the isozymes by relative speed
    speed_enzyme(model, isozyme) =
        max(isozyme.kcat_forward, isozyme.kcat_backward) / sum(
            count * gene_product_molar_mass(model, gid) for
            (gid, count) in isozyme.gene_product_stoichiometry
        )

    # helper function to return the fastest isozyme or nothing
    ris_(model, rid) = begin
        isozymes = reaction_isozymes(model, rid)
        isnothing(isozymes) && return nothing
        argmax(isozyme -> speed_enzyme(model, isozyme), isozymes)
    end

    columns = Vector{Wrappers.Internal.SimplifiedEnzymeConstrainedColumn}()

    (lbs, ubs) = bounds(model)
    rids = variables(model)

    for i = 1:n_variables(model)

        isozyme = ris_(model, rids[i])

        if isnothing(isozyme)
            # non-enzymatic reaction (or a totally ignored one)
            push!(
                columns,
                Wrappers.Internal.SimplifiedEnzymeConstrainedColumn(
                    i,
                    0,
                    lbs[i],
                    ubs[i],
                    0,
                ),
            )
            continue
        end

        mw = sum(
            gene_product_molar_mass(model, gid) * ps for
            (gid, ps) in isozyme.gene_product_stoichiometry
        )

        if min(lbs[i], ubs[i]) < 0 && isozyme.kcat_backward > constants.tolerance
            # reaction can run in reverse
            push!(
                columns,
                Wrappers.Internal.SimplifiedEnzymeConstrainedColumn(
                    i,
                    -1,
                    max(-ubs[i], 0),
                    -lbs[i],
                    mw / isozyme.kcat_backward,
                ),
            )
        end

        if max(lbs[i], ubs[i]) > 0 && isozyme.kcat_forward > constants.tolerance
            # reaction can run forward
            push!(
                columns,
                Wrappers.Internal.SimplifiedEnzymeConstrainedColumn(
                    i,
                    1,
                    max(lbs[i], 0),
                    ubs[i],
                    mw / isozyme.kcat_forward,
                ),
            )
        end
    end

    return SimplifiedEnzymeConstrainedModel(columns, total_gene_product_mass_bound, model)
end
