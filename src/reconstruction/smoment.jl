
"""
$(TYPEDSIGNATURES)

Construct a model with a structure given by sMOMENT algorithm; returns a
[`SMomentModel`](@ref) (see the documentation for details).

# Arguments

- `reaction_isozyme` parameter is a function that returns a single
  [`Isozyme`](@ref) for each reaction, or `nothing` if the reaction is not
  enzymatic. If the reaction has multiple isozymes, use
  [`smoment_isozyme_speed`](@ref) to select the fastest one, as recommended by
  the sMOMENT paper.
- `gene_product_molar_mass` parameter is a function that returns a molar mass
  of each gene product as specified by sMOMENT.
- `total_enzyme_capacity` is the maximum "enzyme capacity" in the model.

Alternatively, all function arguments also accept dictionaries that are used to
provide the same data lookup.
"""
function make_smoment_model(
    model::AbstractMetabolicModel;
    reaction_isozyme::Union{Function,Dict{String,Isozyme}},
    gene_product_molar_mass::Union{Function,Dict{String,Float64}},
    total_enzyme_capacity::Float64,
)
    ris_ =
        reaction_isozyme isa Function ? reaction_isozyme :
        (rid -> get(reaction_isozyme, rid, nothing))
    gpmm_ =
        gene_product_molar_mass isa Function ? gene_product_molar_mass :
        (gid -> gene_product_molar_mass[gid])

    columns = Vector{Types._SMomentColumn}()

    (lbs, ubs) = bounds(model)
    rids = reactions(model)

    for i = 1:n_reactions(model)
        isozyme = ris_(rids[i])
        if isnothing(isozyme)
            # non-enzymatic reaction (or a totally ignored one)
            push!(columns, Types._SMomentColumn(i, 0, lbs[i], ubs[i], 0))
            continue
        end

        mw = sum(gpmm_(gid) * ps for (gid, ps) in isozyme.gene_product_count)

        if min(lbs[i], ubs[i]) < 0 && isozyme.kcat_reverse > constants.tolerance
            # reaction can run in reverse
            push!(
                columns,
                Types._SMomentColumn(
                    i,
                    -1,
                    max(-ubs[i], 0),
                    -lbs[i],
                    mw / isozyme.kcat_reverse,
                ),
            )
        end

        if max(lbs[i], ubs[i]) > 0 && isozyme.kcat_forward > constants.tolerance
            # reaction can run forward
            push!(
                columns,
                Types._SMomentColumn(
                    i,
                    1,
                    max(lbs[i], 0),
                    ubs[i],
                    mw / isozyme.kcat_forward,
                ),
            )
        end
    end

    return SMomentModel(columns, total_enzyme_capacity, model)
end
