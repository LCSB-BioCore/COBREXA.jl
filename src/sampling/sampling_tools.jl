"""
    get_warmup_points(cbm; random_objective=false, numstop=1e10)

Generate warmup points for all the reactions in the model that
are not fixed. Assumes you feed in a JuMP model that is already
constrained i.e. the constraints are already applied to `cbm`.

If you do not want to generate all possible warmup points 
(2*n where n is the number of non-fixed variables), then reduce `numstop`
to the desired number of warmup points. This number of warmup points will
be randomly drawn from the non-fixed fluxes.

By default each flux will be minimized and maximized, however, random objectives
can be used by setting `random_objective` true.
"""
function get_warmup_points(cbm; random_objective = false, numstop = 1e10)
    v = cbm[:x]
    ubs = cbm[:ubs]
    lbs = cbm[:lbs]
    # determine which rxns should be max/min-ized (non fixed fluxes)
    fixed_rxns = Int64[]
    for i in eachindex(v)
        ub_val = normalized_rhs(ubs[i])
        lb_val = normalized_rhs(lbs[i])
        if ub_val >= 0 && ub_val == -lb_val
            push!(fixed_rxns, i)
        elseif ub_val < 0 && -ub_val == lb_val
            push!(fixed_rxns, i)
        end
    end

    # determine number of warmup points (if numstops is less than the number of non-fixed reactions.)
    var_rxn_inds = shuffle!(filter(x -> !(x in fixed_rxns), 1:length(v))) # shuffle the optimization points
    NN = numstop > length(var_rxn_inds) ? length(var_rxn_inds) : numstop
    wpoints = zeros(length(v), 2 * NN)

    # generate warmup points.
    for (i, ii) in zip(1:length(var_rxn_inds), 1:2:(2*length(var_rxn_inds)))
        i > NN && break

        if random_objective
            @objective(cbm, Max, sum(rand() * v[iii] for iii in var_rxn_inds))
        else
            @objective(cbm, Max, v[var_rxn_inds[i]])
        end

        optimize!(cbm)
        for j = 1:size(wpoints, 1)
            wpoints[j, ii] = value(v[j])
        end

        if random_objective
            @objective(cbm, Min, sum(rand() * v[iii] for iii in var_rxn_inds))
        else
            @objective(cbm, Min, v[var_rxn_inds[i]])
        end

        optimize!(cbm)
        for j = 1:size(wpoints, 1)
            wpoints[j, ii+1] = value(v[j])
        end
    end

    return wpoints
end

"""
    get_bound_vectors(opt_model)

Return Float64 vectors of the lower and upper bounds of `opt_model` constraints, 
which is a JuMP model.
"""
function get_bound_vectors(opt_model)
    lbconref = opt_model[:lbs]
    ubconref = opt_model[:ubs]
    lbs = zeros(length(lbconref))
    for i in eachindex(lbs)
        lbval = normalized_rhs(lbconref[i])
        if lbval > 0
            lbs[i] = -abs(lbval)
        else
            lbs[i] = abs(lbval)
        end
    end
    ubs = [normalized_rhs(ubconref[i]) for i in eachindex(ubconref)]

    return lbs, ubs
end

"""
    hit_and_run(N::Int64, cbm; keepevery=100, samplesize=1000, random_objective=false)

Perform basic hit and run sampling for `N` iterations on `cbm`.
Note that `cbm` is a JuMP model that contains whatever solver will be used to find the warmup points.
Also note, `cbm` should already be fully constrained by however you desire.
Every `keepevery` iteration is logged as a sample, where the sample size matrix has `samplesize` columns.
Warm up points are generated in a flux variability sense, unless `random_objective` is true,
in which case a randomly weighted objective is used 2*number of reactions to define the warmup points.

Note that N needs to be >> samplesize.
Sample size is the size of the samples kept in memory.
The larger samplesize is the better the approximation becomes, but the more memory the sampler requires.

See also: [`achr`](@ref)
"""
function hit_and_run(
    N::Int64,
    cbm;
    keepevery = 100,
    samplesize = 1000,
    random_objective = false,
)
    lbs, ubs = get_bound_vectors(cbm) # get actual ub and lb constraints, can't use model function because the user may have changed them in the function arguments

    wpoints = get_warmup_points(cbm; random_objective = random_objective)

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

"""
    test_samples(samples::Array{Float64, 2}, model::StandardModel, ubs, lbs)

Test if samples pass tests: mass balances and constraints are satisfied..
"""
function test_samples(samples::Array{Float64,2}, mass_balance, balance, lbs, ubs)
    violations = Int64[]
    tol = 1e-6
    for i = 1:size(samples, 2)
        if sum(abs, mass_balance * samples[:, i] - balance) < tol
            equality = true
        else
            equality = false
        end
        inequality =
            all(abs.(lbs .- samples[:, i]) .<= tol) .==
            all(abs.(samples[:, i] .- ubs) .<= tol)
        if !all([equality, inequality])
            push!(violations, i)
        end
    end

    return violations
end

# """
#     achr(N::Int64, model::StandardModel, optimizer;constraints=Dict{String, Tuple{Float64,Float64}}(), keepevery=100, samplesize=1000, solver_attributes=Dict{Any, Any}(), random_objective=false)

# Perform artificially centered hit and run.
# Uses the same arguments as the `hit_and_run` sampler.
# Needs work, for long iterations it becomes unstable (violates bounds).

# See also: [`hit_and_run`](@ref)
# """
# function achr(
#     N::Int64,
#     model::StandardModel,
#     optimizer;
#     constraints = Dict{String,Tuple{Float64,Float64}}(),
#     keepevery = 100,
#     samplesize = 1000,
#     solver_attributes = Dict{Any,Any}(),
#     random_objective = false,
# )
#     cbm, v, _, ubcons, lbcons = build_cbm(model)

#     set_optimizer(cbm, optimizer) # choose optimizer
#     if !isempty(solver_attributes) # set other attributes
#         for (k, v) in solver_attributes
#             set_optimizer_attribute(cbm, k, v)
#         end
#     end

#     # set additional constraints
#     for (rxnid, con) in constraints
#         ind = model.reactions[findfirst(model.reactions, rxnid)]
#         set_bound(ind, ubcons, lbcons; ub = con[1], lb = con[2])
#     end

#     ubs, lbs = get_bound_vectors(ubcons, lbcons) # get actual ub and lb constraints

#     wpoints = get_warmup_points(cbm, v, ubcons, lbcons, random_objective = random_objective)

#     nwpts = size(wpoints, 2) # number of warmup points generated
#     samples = zeros(size(wpoints, 1), samplesize) # sample storage
#     current_point = zeros(size(wpoints, 1))
#     current_point .= wpoints[:, rand(1:nwpts)] # pick random initial point

#     shat = mean(wpoints, dims = 2)[:] # mean point

#     δdirtol = 1e-6 # too small directions get ignored ≈ 0 (solver precision issue)
#     sample_num = 0
#     samplelength = 0
#     updatesamplesizelength = true
#     for n = 1:N
#         if updatesamplesizelength # switch to samples
#             direction_point = (@view wpoints[:, rand(1:nwpts)]) - (@view current_point[:]) # use warmup points to find direction in warmup phase
#         else
#             direction_point = (@view samples[:, rand(1:(samplelength))]) - (@view shat[:]) # after warmup phase, only find directions in sampled space
#         end

#         λmax = 1e10
#         λmin = -1e10
#         for i in eachindex(lbs)
#             δlower = lbs[i] - current_point[i]
#             δupper = ubs[i] - current_point[i]
#             if direction_point[i] < -δdirtol
#                 lower = δupper / direction_point[i]
#                 upper = δlower / direction_point[i]
#             elseif direction_point[i] > δdirtol
#                 lower = δlower / direction_point[i]
#                 upper = δupper / direction_point[i]
#             else
#                 lower = -1e10
#                 upper = 1e10
#             end
#             lower > λmin && (λmin = lower) # max min step size that satisfies all bounds
#             upper < λmax && (λmax = upper) # min max step size that satisfies all bounds
#         end

#         if λmax <= λmin || λmin == -1e10 || λmax == 1e10 # this sometimes can happen
#             @warn "Infeasible direction at iteration $(n)..."
#             continue
#         end

#         λ = rand() * (λmax - λmin) + λmin # random step size
#         current_point .= current_point .+ λ .* direction_point # will be feasible
#         shat .= (n .* shat .+ current_point) ./ (n + 1)

#         if n % keepevery == 0
#             sample_num += 1
#             samples[:, sample_num] .= current_point
#             if sample_num >= samplesize
#                 updatesamplesizelength = false
#                 sample_num = 0 # reset, start replacing the older samples
#             end
#             updatesamplesizelength && (samplelength += 1) # lags sample_num because the latter is a flag as well
#         end

#     end

#     violation_inds = test_samples(samples, model, ubs, lbs)
#     if !isempty(violation_inds)
#         @warn "Samples: $violation_inds do not satisfy the problem constraints."
#     end

#     return samples
# end
