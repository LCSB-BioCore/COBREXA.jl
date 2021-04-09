"""
    hit_and_run(N::Int, opt_model; keepevery=100, samplesize=1000, random_objective=false)

Perform basic hit and run sampling for `N` iterations on `opt_model`, where `opt_model` is a constraint based model.
Note that `opt_model` is a JuMP model that contains whatever solver will be used to find the warmup points.
Also note, `opt_model` should already be fully constrained by however you desire.
Every `keepevery` iteration is logged as a sample, where the sample size matrix has `samplesize` columns.
Warm up points are generated in a flux variability sense, unless `random_objective` is true,
in which case a randomly weighted objective is used 2*number of reactions to define the warmup points.

Note that N needs to be >> samplesize.
Sample size is the size of the samples kept in memory.
The larger samplesize is the better the approximation becomes, but the more memory the sampler requires.

# Example
```
using COBREXA
using JuMP
using Tulip

model = read_model("e_coli_core.json")
biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
glucose = findfirst(model.reactions, "EX_glc__D_e")

opt_model = flux_balance_analysis(model, Tulip.Optimizer; 
    modifications=[modify_objective(biomass), 
    modify_constraint(glucose, -12, -12), 
    modify_solver_attribute("IPM_IterationsLimit", 500)])

biomass_index = model[biomass]
λ = JuMP.value(opt_model[:x][biomass_index])
modify_constraint(biomass, 0.99*λ, λ)(model, opt_model)

samples = hit_and_run(100_000, opt_model; keepevery=10, samplesize=5000)
```

See also: [`achr`](@ref)
"""
function hit_and_run(
    N::Int,
    opt_model;
    keepevery = 100,
    samplesize = 1000,
    random_objective = false,
)

    lbs, ubs = get_bound_vectors(opt_model) # get actual ub and lb constraints, can't use model function because the user may have changed them in the function arguments

    wpoints = get_warmup_points(opt_model; random_objective = random_objective)

    nwpts = size(wpoints, 2) # number of warmup points generated
    samples = zeros(size(wpoints, 1), samplesize) # sample storage
    current_point = zeros(size(wpoints, 1))
    current_point .= wpoints[:, rand(1:nwpts)] # pick random initial point

    δdirtol = 1e-6 # too small directions get ignored ≈ 0 (solver precision issue)
    sample_num = 0
    samplelength = 0
    updatesamplesizelength = true
    for n = 1:N

        # direction = random point - current point
        if updatesamplesizelength
            direction_point = (@view wpoints[:, rand(1:nwpts)]) - (@view current_point[:]) # use warmup points to find direction in warmup phase
        else
            direction_point =
                (@view samples[:, rand(1:(samplelength))]) - (@view current_point[:]) # after warmup phase, only find directions in sampled space
        end

        λmax = 1e10
        λmin = -1e10
        for i in eachindex(lbs)
            δlower = lbs[i] - current_point[i]
            δupper = ubs[i] - current_point[i]
            # only consider the step size bound if the direction of travel is non-negligible
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

        if n % keepevery == 0
            sample_num += 1
            samples[:, sample_num] .= current_point
            if sample_num >= samplesize
                updatesamplesizelength = false # once the entire memory vector filled, stop using warm up points
                sample_num = 0 # reset, start replacing the older samples
            end
            updatesamplesizelength && (samplelength += 1) # lags sample_num because the latter is a flag as well
        end

    end

    return samples
end
