
"""
    envelope_lattice(model::MetabolicModel, rids::Vector{String}; kwargs...)

Version of [`envelope_lattice`](@ref) that works on string reaction IDs instead
of integer indexes.
"""
envelope_lattice(model::MetabolicModel, rids::Vector{String}; kwargs...) =
    envelope_lattice(model, Vector{Int}(indexin(rids, reactions(model))); kwargs...)

"""
    envelope_lattice(
        model::MetabolicModel,
        ridxs::Vector{Int};
        samples = 10,
        ranges = collect(zip(bounds(model)...))[ridxs],
        reaction_samples = fill(samples, length(ridxs)),
    )

Create a lattice (list of "tick" vectors) for reactions at indexes `ridxs` in a
model. Arguments `samples`, `ranges`, and `reaction_samples` may be optionally
specified to customize the latice creation process.
"""
envelope_lattice(
    model::MetabolicModel,
    ridxs::Vector{Int};
    samples = 10,
    ranges = collect(zip(bounds(model)...))[ridxs],
    reaction_samples = fill(samples, length(ridxs)),
) = (
    lb .+ (ub - lb) .* ((1:s) .- 1) ./ max(s - 1, 1) for
    (s, (lb, ub)) in zip(reaction_samples, ranges)
)

"""
    objective_envelope(model::MetabolicModel, rids::Vector{String}, args...; kwargs...)

Versioin of [`objective_envelope`](@ref) that works on string reaction IDs
instead of integer indexes.
"""
objective_envelope(model::MetabolicModel, rids::Vector{String}, args...; kwargs...) =
    objective_envelope(
        model,
        Vector{Int}(indexin(rids, reactions(model))),
        args...;
        kwargs...,
    )

"""
    objective_envelope(
        model::MetabolicModel,
        ridxs::Vector{Int},
        optimizer;
        lattice = envelope_lattice(model, ridxs),
        kwargs...,
    )

Compute an array of objective values for the `model` for rates of reactions
specified `ridxs` fixed to a regular range of values between their respective
lower and upper bounds.

This can be used to compute a "production envelope" of a metabolite; but
generalizes to any specifiable objective and to multiple dimensions of the
examined space. To retrieve a production envelope of any metabolite, set the
objective coefficient vector of the `model` to a vector that contains a single
`1` for the exchange reaction that "outputs" this metabolite, and run
[`objective_envelope`](@ref) with the exchange reaction of the "parameter"
metabolite specified in `ridxs`.

Returns a named tuple that contains `lattice` with reference values of the
metabolites, and an N-dimensional array `values` with the computed objective
values, where N is the number of specified reactions.  Because of the
increasing dimensionality, the computation gets very voluminous with increasing
length of `ridxs`.

`kwargs` are internally forwarded to [`screen_optmodel_modifications`](@ref).
"""
objective_envelope(
    model::MetabolicModel,
    ridxs::Vector{Int},
    optimizer;
    lattice = envelope_lattice(model, ridxs),
    kwargs...,
) = (
    lattice = collect.(lattice),
    values = screen_optmodel_modifications(
        model,
        optimizer;
        modifications = collect(
            [(_, optmodel) -> begin
                    for (i, ridx) in enumerate(ridxs)
                        set_normalized_rhs(optmodel[:lbs][ridx], fluxes[i])
                        set_normalized_rhs(optmodel[:ubs][ridx], fluxes[i])
                    end
                end] for fluxes in Iterators.product(lattice...)
        ),
        analysis = screen_optimize_objective,
        kwargs...,
    ),
)
