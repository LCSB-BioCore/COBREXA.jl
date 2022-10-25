"""
$(TYPEDSIGNATURES)

Run a hit-and-run style sampling that starts from `warmup_points` and uses
their affine combinations for generating the run directions to sample the space
delimited by `lbs` and `ubs`.  The reaction rate vectors in `warmup_points`
should be organized in columns, i.e. `warmup_points[:,1]` is the first set of
reaction rates.

There are total `chains` of hit-and-run runs, each on a batch of
`size(warmup_points, 2)` points. The runs are scheduled on `workers`, for good
load balancing `chains` should be ideally much greater than `length(workers)`.

Each run continues for `maximum(sample_iters)` iterations; the numbers in
`sample_iters` represent the iterations at which the whole "current" batch of
points is collected for output. For example, `sample_iters=[1,4,5]` causes the
process run for 5 iterations, returning the sample batch that was produced by
1st, 4th and last (5th) iteration.

Returns a matrix of sampled reaction rates (in columns), with all collected
samples horizontally concatenated. The total number of samples (columns) will
be `size(warmup_points,2) * chains * length(sample_iters)`.

# Example
```
warmup_points = warmup_from_variability(model, GLPK.Optimizer)
samples = affine_hit_and_run(model, warmup_points, sample_iters = 101:105)

# convert the result to flux (for models where the distinction matters):
fluxes = reaction_flux(model)' * samples
```
"""
function affine_hit_and_run(
    m::AbstractMetabolicModel,
    warmup_points::Matrix{Float64};
    sample_iters = 100 .* (1:5),
    workers = [myid()],
    chains = length(workers),
    seed = rand(Int),
)
    @assert size(warmup_points, 1) == n_reactions(m)

    lbs, ubs = bounds(m)
    C = coupling(m)
    cl, cu = coupling_bounds(m)
    if isnothing(C)
        C = zeros(0, n_reactions(m))
        cl = zeros(0)
        cu = zeros(0)
    end
    save_at.(workers, :cobrexa_hit_and_run_data, Ref((warmup_points, lbs, ubs, C, cl, cu)))

    # sample all chains
    samples = hcat(
        pmap(
            chain -> _affine_hit_and_run_chain(
                (@remote cobrexa_hit_and_run_data)...,
                sample_iters,
                seed + chain,
            ),
            CachingPool(workers),
            1:chains,
        )...,
    )

    # remove warmup points from workers
    map(fetch, remove_from.(workers, :cobrexa_hit_and_run_data))

    return samples
end

"""
$(TYPEDSIGNATURES)

Internal helper function for computing a single affine hit-and-run chain.
"""
function _affine_hit_and_run_chain(warmup, lbs, ubs, C, cl, cu, iters, seed)

    rng = StableRNG(seed % UInt)
    points = copy(warmup)
    d, n_points = size(points)
    n_couplings = size(C, 1)
    result = Matrix{Float64}(undef, d, n_points * length(iters))

    # helper for reducing the available run range
    function update_range(range, pos, dir, lb, ub)
        dl = lb - pos
        du = ub - pos
        lower, upper =
            dir < -constants.tolerance ? (du, dl) ./ dir :
            dir > constants.tolerance ? (dl, du) ./ dir : (-Inf, Inf)
        return (max(range[1], lower), min(range[2], upper))
    end

    iter = 0

    for (iter_idx, iter_target) in enumerate(iters)

        while iter < iter_target
            iter += 1

            new_points = copy(points)

            for i = 1:n_points

                mix = rand(rng, n_points) .+ constants.tolerance
                dir = points * (mix ./ sum(mix)) - points[:, i]

                # iteratively collect the maximum and minimum possible multiple
                # of `dir` added to the current point
                run_range = (-Inf, Inf)
                for j = 1:d
                    run_range =
                        update_range(run_range, points[j, i], dir[j], lbs[j], ubs[j])
                end

                # do the same for coupling
                dc = C * dir
                pc = C * points[:, i]
                for j = 1:n_couplings
                    run_range = update_range(run_range, pc[j], dc[j], cl[j], cu[j])
                end

                # generate a point in the viable run range and update it
                lambda = run_range[1] + rand(rng) * (run_range[2] - run_range[1])
                isfinite(lambda) || continue # avoid divergence
                new_points[:, i] = points[:, i] .+ lambda .* dir

                # TODO normally, here we would check if sum(S*new_point) is still
                # lower than the tolerance, but we shall trust the computer
                # instead.
            end

            points = new_points
        end

        result[:, n_points*(iter_idx-1)+1:iter_idx*n_points] .= points
    end

    result
end
