"""
$(TYPEDSIGNATURES)

Wrap a `model` with a structure given by sMOMENT algorithm. Returns a
[`SimplifiedEnzymeConstrainedModel`](@ref) (see the documentation for details).
The sMOMENT algorithm only uses one [`Isozyme`](@ref) per reaction. If multiple
isozymes are present in `model`, the "fastest" isozyme will be used. This is
determined based on maximum kcat (forward or backward) divided by mass of the
isozyme. Multiple enzyme capacity constraint can be placed on the model using
the keyword arguments.

Parameters `reaction_mass_groups` and `reaction_mass_group_bounds` specify
groups of reactions (by their IDs), and their respective total mass limit.
Reactions that are not listed in any reaction mass group are ignored (likewise
if they don't have isozymes).

For simplicity, specifying the `total_reaction_mass_bound` argument overrides
the above arguments by internally specifying a single group on all reactions
that acts like a maximum "enzyme capacity" in the model.

# Example
```
ecmodel = make_simplified_enzyme_constrained_model(
    model;
    reaction_mass_groups = Dict(
        "membrane" => ["r1", "r2"],
        "total" => ["r1", "r2", "r3"],
    ),
    reaction_mass_group_bounds = Dict(
        "membrane" => 0.2,
        "total" => 0.5,
    ),
)

ecmodel2 = make_simplified_enzyme_constrained_model(
    model;
    total_reaction_mass_bound = 0.5
)
```
"""
function make_simplified_enzyme_constrained_model(
    model::AbstractMetabolicModel;
    total_reaction_mass_bound::Maybe{Float64} = nothing,
    reaction_mass_groups::Maybe{Dict{String,Vector{String}}} = nothing,
    reaction_mass_group_bounds::Maybe{Dict{String,Float64}} = nothing,
)
    all(
        !isnothing,
        [reaction_mass_groups, reaction_mass_group_bounds, total_reaction_mass_bound],
    ) && throw(ArgumentError("Too many arguments specified!"))

    # fix kwarg inputs
    if !isnothing(total_reaction_mass_bound)
        reaction_mass_groups = Dict("uncategorized" => variables(model)) # TODO should be reactions
        reaction_mass_group_bounds = Dict("uncategorized" => total_reaction_mass_bound)
    end
    isnothing(reaction_mass_groups) &&
        throw(ArgumentError("missing reaction mass group specification"))
    isnothing(reaction_mass_group_bounds) &&
        throw(ArgumentError("missing reaction mass group bounds"))

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

    bound_ids = keys(reaction_mass_group_bounds)
    total_reaction_mass_bounds = collect(values(reaction_mass_group_bounds))

    (lbs, ubs) = bounds(model) # TODO need a reaction_bounds accessor for full generality
    rids = variables(model) # TODO needs to be reactions

    for i = 1:n_variables(model) # TODO this should be reactions

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
                    Int64[],
                ),
            )
            continue
        end

        mw = sum(
            gene_product_molar_mass(model, gid) * ps for
            (gid, ps) in isozyme.gene_product_stoichiometry
        )

        bidxs = [
            idx for
            (idx, bid) in enumerate(bound_ids) if rids[i] in reaction_mass_groups[bid]
        ]

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
                    bidxs,
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
                    bidxs,
                ),
            )
        end
    end

    return SimplifiedEnzymeConstrainedModel(columns, total_reaction_mass_bounds, model)
end
