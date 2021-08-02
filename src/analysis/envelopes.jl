
envelope_lattice(m::MetabolicModel, rids::Vector{String}; kwargs...) =
    envelope_lattice(m, Vector{Int}(indexin(rids, reactions(m))); kwargs...)

envelope_lattice(
    m::MetabolicModel,
    ridxs::Vector{Int};
    samples = 10,
    ranges = collect(zip(bounds(m)...))[ridxs],
    reaction_samples = fill(samples, length(ridxs))) =
    (
        lb .+ (ub - lb) .* ((1:s) .- 1) ./ max(s - 1, 1) for
        (s, (lb, ub)) in zip(reaction_samples, ranges)
    )

objective_envelope(m::MetabolicModel, rids::Vector{String}, args...; kwargs...) =
    objective_envelope(m, Vector{Int}(indexin(rids, reactions(m))), args...; kwargs...)

objective_envelope(
    m::MetabolicModel,
    ridxs::Vector{Int},
    optimizer;
    lattice = envelope_lattice(m, ridxs),
    kwargs...,
) = (
    lattice = collect.(lattice),
    values = screen_optmodel_modifications(
        m,
        optimizer,
        modifications = collect(
            [(_, optmodel) -> begin
                    for (i, ridx) in enumerate(ridxs)
                        set_normalized_rhs(optmodel[:lbs][ridx], bounds[i])
                        set_normalized_rhs(optmodel[:ubs][ridx], bounds[i])
                    end
                end] for bounds in Iterators.product(lattice...)
        ),
        analysis = screen_optimize_objective,
    ),
)
