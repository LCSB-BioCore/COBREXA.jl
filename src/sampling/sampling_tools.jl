"""
    get_warmup_points(cbmodel, v, mb, lbs, ubs; random_objective=false, numstop=1e10)

Generate warmup points for all the reactions on the model that 
are not fixed. Assumes you feed in a JuMP model that is already
constrained i.e. the constrains are already applied into cbmodel.
Note, extra constraints applied to ubs and lbs will have no effect.

numstop = 2*number of warmup points - to reduce the time this takes
"""
function get_warmup_points(cbmodel, v, lbs, ubs; random_objective = false, numstop = 1e10)
    # determine which rxns should be max/min-ized
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

    # determine number of warmup points
    var_rxn_inds = shuffle!(filter(x -> !(x in fixed_rxns), 1:length(v))) # shuffle the optimization points
    NN = numstop > length(var_rxn_inds) ? length(var_rxn_inds) : numstop
    wpoints = zeros(length(v), 2 * NN)

    for (i, ii) in zip(1:length(var_rxn_inds), 1:2:(2*length(var_rxn_inds)))
        i > NN && break

        if random_objective
            @objective(cbmodel, Max, sum(rand() * v[iii] for iii in var_rxn_inds))
        else
            @objective(cbmodel, Max, v[var_rxn_inds[i]])
        end

        optimize!(cbmodel)
        for j = 1:size(wpoints, 1)
            wpoints[j, ii] = value(v[j])
        end

        if random_objective
            @objective(cbmodel, Min, sum(rand() * v[iii] for iii in var_rxn_inds))
        else
            @objective(cbmodel, Min, v[var_rxn_inds[i]])
        end

        optimize!(cbmodel)
        for j = 1:size(wpoints, 1)
            wpoints[j, ii+1] = value(v[j])
        end
    end

    return wpoints
end

"""
    get_bound_vectors(ubconref, lbconref)

Return Float64 vectors of the lower and upper bounds of the JuMP constraint refs.
"""
function get_bound_vectors(lbconref, ubconref)
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
    hit_and_run(N::Int64, model::CobraModel, optimizer; constraints=Dict{String, Tuple{Float64,Float64}}(), keepevery=100, samplesize=1000, solver_attributes=Dict{Any, Any}(), random_objective=false)

Perform basic hit and run sampling for `N` iterations using `model` with `optimizer` from `JuMP`. 
Additional constraints supplied by `constraints` as a dictionary of reaction `id`s mapped to a tuple of `(lb, ub)` of fluxes.
Every `keepevery` iteration is logged as a sample, where the sample size matrix has `samplesize` columns.
Solver specific settings can be set using `solver_attributes`.
Warm up points are generated in a flux variability sense, unless `random_objective` is true, in which case a randomly weighted objective is used 2*number of reactions to define the warmup points.

Note that N needs to be >> samplesize. 
Sample size is the size of the samples kept in memory. 
The larger samplesize is the better the approximation becomes, but the more memory the sampler requires.

See also: [`achr`](@ref)
"""
function hit_and_run(
    N::Int64,
    model::CobraModel,
    optimizer;
    constraints = Dict{String,Tuple{Float64,Float64}}(),
    keepevery = 100,
    samplesize = 1000,
    solver_attributes = Dict{Any,Any}(),
    random_objective = false,
    sense = MOI.MAX_SENSE,
)
    # get core optimization problem
    cbmodel, v, mb, lbcons, ubcons =
        make_optimization_model(model, optimizer, sense = sense)

    if !isempty(solver_attributes) # set other attributes
        for (k, v) in solver_attributes
            set_optimizer_attribute(cbmodel, k, v)
        end
    end

    # set additional constraints
    for (rxnid, con) in constraints
        ind = model.reactions[findfirst(model.reactions, rxnid)]
        set_bound(ind, lbcons, ubcons; lb = con[1], ub = con[2])
    end

    lbs, ubs = get_bound_vectors(lbcons, ubcons) # get actual ub and lb constraints, can't use model function because the user may have changed them in the function arguments

    wpoints =
        get_warmup_points(cbmodel, v, lbcons, ubcons, random_objective = random_objective)

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

    violation_inds = test_samples(samples, model, ubs, lbs)
    if !isempty(violation_inds)
        @warn "Samples: $violation_inds do not satisfy the problem constraints."
    end

    return samples
end

"""
    test_samples(samples::Array{Float64, 2}, model::CobraModel, ubs, lbs)

Test if samples pass tests: mass balances and constraints are satisfied..
"""
function test_samples(samples::Array{Float64,2}, model::CobraModel, ubs, lbs)
    S = stoichiometry(model) # assume S has not been modified from model AND b is zero
    violations = Int64[]
    tol = 1e-6
    for i = 1:size(samples, 2)
        if sum(abs, S * samples[:, i]) < tol
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

"""
    achr(N::Int64, model::CobraModel, optimizer;constraints=Dict{String, Tuple{Float64,Float64}}(), keepevery=100, samplesize=1000, solver_attributes=Dict{Any, Any}(), random_objective=false)

Perform artificially centered hit and run.
Uses the same arguments as the `hit_and_run` sampler.
Needs work, for long iterations it becomes unstable (violates bounds).

See also: [`hit_and_run`](@ref)
"""
function achr(
    N::Int64,
    model::CobraModel,
    optimizer;
    constraints = Dict{String,Tuple{Float64,Float64}}(),
    keepevery = 100,
    samplesize = 1000,
    solver_attributes = Dict{Any,Any}(),
    random_objective = false,
)
    cbmodel, v, _, ubcons, lbcons = build_cbm(model)

    set_optimizer(cbmodel, optimizer) # choose optimizer
    if !isempty(solver_attributes) # set other attributes
        for (k, v) in solver_attributes
            set_optimizer_attribute(cbmodel, k, v)
        end
    end

    # set additional constraints
    for (rxnid, con) in constraints
        ind = model.reactions[findfirst(model.reactions, rxnid)]
        set_bound(ind, ubcons, lbcons; ub = con[1], lb = con[2])
    end

    ubs, lbs = get_bound_vectors(ubcons, lbcons) # get actual ub and lb constraints

    wpoints =
        get_warmup_points(cbmodel, v, ubcons, lbcons, random_objective = random_objective)

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
