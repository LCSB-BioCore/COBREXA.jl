"""
    function affine_hit_and_run(
        warmup_points::Matrix{Float64},
        lbs::Vector{Float64},
        ubs::Vector{Float64};
        sample_iters = 100 .* (1:5),
        workers = [myid()],
        chains = length(workers),
        seed = rand(Int),
    )

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
using COBREXA
using Tulip

model = load_model(StandardModel, model_path)

warmup, lbs, ubs = warmup_from_variability(model, Tulip.Optimizer, 100)
samples = affine_hit_and_run(warmup, lbs, ubs, sample_iters = 1:3)

# convert the result to flux (for models where the distinction matters):
fluxes = reaction_flux(model)' * samples
```
"""
function affine_hit_and_run(
    warmup_points::Matrix{Float64},
    lbs::Vector{Float64},
    ubs::Vector{Float64};
    sample_iters = 100 .* (1:5),
    workers = [myid()],
    chains = length(workers),
    seed = rand(Int),
)

    # distribute starting data to workers
    save_at.(workers, :cobrexa_hit_and_run_data, Ref((warmup_points, lbs, ubs)))

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
    _affine_hit_and_run_chain(warmup, lbs, ubs, iters, seed)

Internal helper function for computing a single affine hit-and-run chain.
"""
function _affine_hit_and_run_chain(warmup, lbs, ubs, iters, seed)

    rng = StableRNG(seed % UInt)
    points = copy(warmup)
    d, n_points = size(points)
    result = Matrix{Float64}(undef, size(points, 1), 0)

    iter = 0

    for iter_target in iters

        while iter < iter_target
            iter += 1

            new_points = copy(points)

            for i = 1:n_points

                mix = rand(rng, n_points) .+ _constants.tolerance
                dir = points * (mix ./ sum(mix)) - points[:, i]

                # iteratively collect the maximum and minimum possible multiple
                # of `dir` added to the current point
                lambda_max = Inf
                lambda_min = -Inf
                for j = 1:d
                    dl = lbs[j] - points[j, i]
                    du = ubs[j] - points[j, i]
                    idir = 1 / dir[j]
                    if dir[j] < -_constants.tolerance
                        lower = du * idir
                        upper = dl * idir
                    elseif dir[j] > _constants.tolerance
                        lower = dl * idir
                        upper = du * idir
                    else
                        lower = -Inf
                        upper = Inf
                    end
                    lambda_min = max(lambda_min, lower)
                    lambda_max = min(lambda_max, upper)
                end

                lambda = lambda_min + rand(rng) * (lambda_max - lambda_min)
                !isfinite(lambda) && continue # avoid divergence
                new_points[:, i] = points[:, i] .+ lambda .* dir

                # TODO normally, here we would check if sum(S*new_point) is still
                # lower than the tolerance, but we shall trust the computer
                # instead.
            end

            points = new_points
        end

        result = hcat(result, points)
    end

    result
end
