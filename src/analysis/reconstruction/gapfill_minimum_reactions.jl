"""
$(TYPEDSIGNATURES)

Find a minimal set of reactions from `universal_reactions` that should be added
to `model` so that the model has a feasible solution with bounds on its
objective function given in `objective_bounds`.  Weights of the added reactions
may be specified in `weights` to prefer adding reactions with lower weights.

Internally, this builds and solves a mixed integer program, following
the method of Reed et al. (Reed, Jennifer L., et al. "Systems approach to
refining genome annotation." *Proceedings of the National Academy of Sciences*
(2006)).

The function returns a solved JuMP optimization model, with the boolean
reaction inclusion indicators in variable vector `y`. Use
[`gapfilled_mask`](@ref) or [`gapfilled_rids`](@ref) to collect the reaction
information in Julia datatypes.

To reduce the uncertainty in the MILP solver (and likely reduce the
complexity), you may put a limit on the size of the added reaction set in
`maximum_new_reactions`.

# Common pitfalls

If [`gapfill_minimum_reactions`](@ref) is supposed to generate any reasonable
output, the input model *MUST NOT* be feasible, otherwise there is "no work to
do" and no reactions are added. Notably, an inactive model (the flux is zero)
is considered to be feasible. If this is the case, [`gapfilled_rids`](@ref)
will return an empty vector (as opposed to `nothing`).

To prevent this, you may need to modify the model to disallow the trivial
solutions (for example by putting a lower bound on reactions that you expect to
be working in the solved model, in a similar manner like how the ATP
maintenance reaction is bounded in E. Coli "core" model). The
`objective_bounds` parameter makes this easier by directly placing a bound on
the objective value of the model, which typically forces the model to be
active.

The `maximum_new_reactions` parameter may have critical impact on performance
in some solvers, because (in a general worst case) there is
`2^maximum_new_reactions` model variants to be examined. Putting a hard limit
on the reaction count may serve as a heuristic that helps the solver not to
waste too much time solving impractically complex subproblems.
"""
function gapfill_minimum_reactions(
    model::AbstractMetabolicModel,
    universal_reactions::Vector{Reaction},
    optimizer;
    objective_bounds = (constants.tolerance, constants.default_reaction_bound),
    maximum_new_reactions = length(universal_reactions),
    weights = fill(1.0, length(universal_reactions)),
    modifications = [],
)
    precache!(model)

    # constraints from universal reactions that can fill gaps
    univs = _universal_stoichiometry(universal_reactions, metabolites(model))

    # add space for additional metabolites and glue with the universal reaction
    # stoichiometry
    extended_stoichiometry = [[
        stoichiometry(model)
        spzeros(length(univs.new_mids), n_reactions(model))
    ] univs.stoichiometry]

    # make the model anew (we can't really use make_optimization_model because
    # we need the balances and several other things completely removed. Could
    # be solved either by parametrizing make_optimization_model or by making a
    # tiny temporary wrapper for this.
    # keep this in sync with src/base/solver.jl, except for adding balances.
    opt_model = Model(optimizer)
    @variable(opt_model, x[1:n_reactions(model)])
    xl, xu = bounds(model)
    @constraint(opt_model, lbs, xl .<= x)
    @constraint(opt_model, ubs, x .<= xu)

    C = coupling(model)
    isempty(C) || begin
        cl, cu = coupling_bounds(model)
        @constraint(opt_model, c_lbs, cl .<= C * x)
        @constraint(opt_model, c_ubs, C * x .<= cu)
    end

    # add the variables for new stuff
    @variable(opt_model, ux[1:length(universal_reactions)]) # fluxes from universal reactions
    @variable(opt_model, y[1:length(universal_reactions)], Bin) # indicators

    # combined metabolite balances
    @constraint(
        opt_model,
        extended_stoichiometry * [x; ux] .==
        [balance(model); zeros(length(univs.new_mids))]
    )

    # objective bounds
    @constraint(opt_model, objective_bounds[1] <= objective(model)' * x)
    @constraint(opt_model, objective_bounds[2] >= objective(model)' * x)

    # flux bounds of universal reactions with indicators
    @constraint(opt_model, ulb, univs.lbs .* y .<= ux)
    @constraint(opt_model, uub, univs.ubs .* y .>= ux)

    # minimize the total number of indicated reactions
    @objective(opt_model, Min, weights' * y)

    # limit the number of indicated reactions
    # (prevents the solver from exploring too far)
    @constraint(opt_model, sum(y) <= maximum_new_reactions)

    # apply all modifications
    for mod in modifications
        mod(model, opt_model)
    end

    optimize!(opt_model)

    return opt_model
end

"""
$(TYPEDSIGNATURES)

Get a `BitVector` of added reactions from the model solved by
[`gapfill_minimum_reactions`](@ref). The bit indexes correspond to the indexes
of `universal_reactions` given to the gapfilling function. In case the model is
not solved, this returns `nothing`.

If this function returns a zero vector (instead of `nothing`), it is very
likely that the original model was already feasible and you may need to
constraint it more. Refer to "pitfalls" section in the documentation of
[`gapfill_minimum_reactions`](@ref) for more details.

# Example

    gapfill_minimum_reactions(myModel, myReactions, Tulip.Optimizer) |> gapfilled_mask
"""
gapfilled_mask(opt_model)::BitVector =
    is_solved(opt_model) ? value.(opt_model[:y]) .> 0 : nothing

"""
$(TYPEDSIGNATURES)

Utility to extract a short vector of IDs of the reactions added by the
gapfilling algorithm. Use with `opt_model` returned from
[`gapfill_minimum_reactions`](@ref).

If this function returns an empty vector (instead of `nothing`), it is very
likely that the original model was already feasible and you may need to
constraint it more. Refer to "pitfalls" section in the documentation of
[`gapfill_minimum_reactions`](@ref) for more details.
"""
gapfilled_rids(opt_model, universal_reactions::Vector{Reaction}) =
    let v = gapfilled_mask(opt_model)
        isnothing(v) ? nothing : [rxn.id for rxn in universal_reactions[v]]
    end

"""
$(TYPEDSIGNATURES)

Overload of [`gapfilled_rids`](@ref) that can be piped easily.

# Example

    gapfill_minimum_reactions(myModel, myReactions, Tulip.Optimizer) |> gapfilled_rids(myReactions)
"""
gapfilled_rids(universal_reactions::Vector{Reaction}) =
    opt_model -> gapfilled_rids(opt_model, universal_reactions)

"""
$(TYPEDSIGNATURES)

A helper function that constructs the stoichiometric matrix of a set of
`universal_reactions`. The order of the metabolites is determined with
`mids`, so that this stoichiometric matrix can be combined with
another one.
"""
function _universal_stoichiometry(urxns::Vector{Reaction}, mids::Vector{String})

    # traversal over all elements in stoichiometry of universal_reactions
    stoiMap(f) = [
        f(ridx, mid, stoi) for (ridx, rxn) in enumerate(urxns) for
        (mid, stoi) in rxn.metabolites
    ]

    # make an index and find new metabolites
    met_id_lookup = Dict(mids .=> eachindex(mids))

    new_mids =
        collect(Set(filter(x -> !haskey(met_id_lookup, x), stoiMap((_, mid, _) -> mid))))
    all_mids = vcat(mids, new_mids)

    # remake the index with all metabolites
    met_id_lookup = Dict(all_mids .=> eachindex(all_mids))

    # build the result
    return (
        stoichiometry = float.(
            sparse(
                stoiMap((_, mid, _) -> met_id_lookup[mid]),
                stoiMap((ridx, _, _) -> ridx),
                stoiMap((_, _, stoi) -> stoi),
                length(all_mids),
                length(urxns),
            ),
        ),
        lbs = [rxn.lower_bound for rxn in urxns],
        ubs = [rxn.upper_bound for rxn in urxns],
        new_mids = new_mids,
    )
end
