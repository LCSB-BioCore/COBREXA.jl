"""
    achr(N::Int64, model::StandardModel, optimizer;constraints=Dict{String, Tuple{Float64,Float64}}(), keepevery=100, samplesize=1000, solver_attributes=Dict{Any, Any}(), random_objective=false)

Perform artificially centered hit and run.
Uses the same arguments as the `hit_and_run` sampler.
Needs work, for long iterations it becomes unstable (violates bounds).

See also: [`hit_and_run`](@ref)
"""
function achr(
    N::Int64,
    model::StandardModel,
    optimizer;
    constraints = Dict{String,Tuple{Float64,Float64}}(),
    keepevery = 100,
    samplesize = 1000,
    solver_attributes = Dict{Any,Any}(),
    random_objective = false,
)
    opt_model, v, _, ubcons, lbcons = build_cbm(model)

    set_optimizer(opt_model, optimizer) # choose optimizer
    if !isempty(solver_attributes) # set other attributes
        for (k, v) in solver_attributes
            set_optimizer_attribute(opt_model, k, v)
        end
    end

    # set additional constraints
    for (rxnid, con) in constraints
        ind = model.reactions[findfirst(model.reactions, rxnid)]
        set_bound(ind, ubcons, lbcons; ub = con[1], lb = con[2])
    end

    ubs, lbs = get_bound_vectors(ubcons, lbcons) # get actual ub and lb constraints

    wpoints =
        get_warmup_points(opt_model, v, ubcons, lbcons, random_objective = random_objective)

    nwpts = size(wpoints, 2) # number of warmup points generated
    samples = zeros(size(wpoints, 1), samplesize) # sample storage
    current_point = zeros(size(wpoints, 1))
    current_point .= wpoints[:, rand(1:nwpts)] # pick random initial point

    shat = mean(wpoints, dims = 2)[:] # mean point

    δdirtol = 1e-6 # too small directions get ignored ≈ 0 (solver precision issue)
    sample_num = 0
    samplelength = 0
    updatesamplesizelength = true
    for n = 1:N
        if updatesamplesizelength # switch to samples
            direction_point = (@view wpoints[:, rand(1:nwpts)]) - (@view current_point[:]) # use warmup points to find direction in warmup phase
        else
            direction_point = (@view samples[:, rand(1:(samplelength))]) - (@view shat[:]) # after warmup phase, only find directions in sampled space
        end

        λmax = 1e10
        λmin = -1e10
        for i in eachindex(lbs)
            δlower = lbs[i] - current_point[i]
            δupper = ubs[i] - current_point[i]
            if direction_point[i] < -δdirtol
                lower = δupper / direction_point[i]
                upper = δlower / direction_point[i]
            elseif direction_point[i] > δdirtol
                lower = δlower / direction_point[i]
                upper = δupper / direction_point[i]
            else
                lower = -1e10
                upper = 1e10
            end
            lower > λmin && (λmin = lower) # max min step size that satisfies all bounds
            upper < λmax && (λmax = upper) # min max step size that satisfies all bounds
        end

        if λmax <= λmin || λmin == -1e10 || λmax == 1e10 # this sometimes can happen
            @warn "Infeasible direction at iteration $(n)..."
            continue
        end

        λ = rand() * (λmax - λmin) + λmin # random step size
        current_point .= current_point .+ λ .* direction_point # will be feasible
        shat .= (n .* shat .+ current_point) ./ (n + 1)

        if n % keepevery == 0
            sample_num += 1
            samples[:, sample_num] .= current_point
            if sample_num >= samplesize
                updatesamplesizelength = false
                sample_num = 0 # reset, start replacing the older samples
            end
            updatesamplesizelength && (samplelength += 1) # lags sample_num because the latter is a flag as well
        end

    end

    violation_inds = test_samples(samples, model, ubs, lbs)
    if !isempty(violation_inds)
        @warn "Samples: $violation_inds do not satisfy the problem constraints."
    end

    return samples
end
