"""
get_warmup_points(cbmodel, v, mb, ubs, lbs; randomobjective=false, numstop=1e10)

Generate warmup points for all the reactions on the model that 
are not fixed. Assumes you feed in a JuMP model that is already
constrained by however you want it to be.

numstop = 2*number of warmup points - to reduce the time this takes
"""
function get_warmup_points(cbmodel, v, ubs, lbs; randomobjective=false, numstop=1e10)
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
    var_rxn_inds =shuffle!(filter(x->!(x in fixed_rxns), 1:length(v))) # shuffle the optimization points
    NN = numstop > length(var_rxn_inds) ? length(var_rxn_inds) : numstop
    wpoints = zeros(length(v), 2*NN)

    for (i, ii) in zip(1:length(var_rxn_inds), 1:2:(2*length(var_rxn_inds)))
        i > NN && break

        if randomobjective
            @objective(cbmodel, Max, sum(rand()*v[iii] for iii in var_rxn_inds))
        else
            @objective(cbmodel, Max, v[var_rxn_inds[i]])
        end

        optimize!(cbmodel)
        for j=1:size(wpoints, 1)
            wpoints[j, ii] = value(v[j])
        end

        if randomobjective
            @objective(cbmodel, Min, sum(rand()*v[iii] for iii in var_rxn_inds))
        else
            @objective(cbmodel, Min, v[var_rxn_inds[i]])
        end

        optimize!(cbmodel)
        for j=1:size(wpoints, 1)
            wpoints[j, ii+1] = value(v[j])
        end
    end

    return wpoints
end

"""
ubs, lbs = get_bound_vectors(ubconref, lbconref)

Return Float64 vectors of the upper and lower bounds of the JuMP
constraint refs.
"""
function get_bound_vectors(ubconref, lbconref)
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
    return ubs, lbs
end

"""
hit_and_run(N::Int64, wpoints::Array{Float64, 2}, model::Model, ubcons, lbcons; keepevery=10, samplesize=1000)

Perform basic hit and run sampling for N iterations. Note that N needs to be >> samplesize. Sample size is the
size of the samples kept in memory. The larger samplesize is the better the approximation becomes, but the more
memory the sampler requires.
"""
function hit_and_run(N::Int64, wpoints::Array{Float64, 2}, ubcons, lbcons; keepevery=100, samplesize=1000)  
    ubs, lbs = get_bound_vectors(ubcons, lbcons) # get bounds from model
    nwpts = size(wpoints, 2) # number of warmup points generated
    samples = zeros(size(wpoints, 1), samplesize) # sample storage
    current_point = zeros(size(wpoints, 1))
    current_point .= wpoints[:, rand(1:nwpts)]

    δdirtol = 1e-6 # too small directions get ignored ≈ 0 (solver precision issue) 
    sample_num = 0
    samplelength = 0
    updatesamplesizelength = true
    for n=1:N

        if updatesamplesizelength
            direction_point = (@view wpoints[:, rand(1:nwpts)]) - (@view current_point[:]) # use warmup points to find direction in warmup phase
        else
            direction_point = (@view samples[:, rand(1:(samplelength))]) - (@view current_point[:]) # after warmup phase, only find directions in sampled space
        end

        λmax = 1e10
        λmin = -1e10 
        for i in eachindex(lbs)
            δlower = lbs[i] - current_point[i]
            δupper = ubs[i] - current_point[i]
            if direction_point[i] < -δdirtol
                lower = δupper/direction_point[i]
                upper = δlower/direction_point[i]
            elseif direction_point[i] > δdirtol
                lower = δlower/direction_point[i]
                upper = δupper/direction_point[i]
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
    
        λ = rand()*(λmax - λmin) + λmin # random step size
        current_point .= current_point .+ λ .* direction_point # will be feasible

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

    return samples
end

"""
indexes_test_failed = test_samples(samples::Array{Float64, 2}, model::Model, ubcons, lbcons)

Test if samples pass tests: mass balances and constraints are satisfied..
"""
function test_samples(samples::Array{Float64, 2}, model::Model, ubcons, lbcons)
    S, _, _, _ = get_core_model(model) # assume S has not been modified from model
    ubs, lbs = get_bound_vectors(ubcons, lbcons)
    violations = Int64[]
    tol = 1e-6
    for i in 1:size(samples, 2)
        if sum(abs, S*samples[:, i]) < tol
            equality = true
        else
            equality = false    
        end
        inequality = all(abs.(lbs .- samples[:, i]) .<= tol) .== all(abs.(samples[:, i] .- ubs) .<= tol) 
        if !all([equality, inequality])
            push!(violations, i)
        end
    end
    
    return violations
end

"""
achr(N::Int64, wpoints::Array{Float64, 2}, ubcons, lbcons; keepevery=100, samplesize=1000)

Needs work, for long iterations it becomes unstable (violates bounds).
"""
function achr(N::Int64, wpoints::Array{Float64, 2}, ubcons, lbcons; keepevery=100, samplesize=1000)  
    ubs, lbs = get_bound_vectors(ubcons, lbcons) # get bounds from model
    nwpts = size(wpoints, 2) # number of warmup points generated
    samples = zeros(size(wpoints, 1), samplesize) # sample storage
    current_point = zeros(size(wpoints, 1))
    current_point .= wpoints[:, rand(1:nwpts)]
    shat = mean(wpoints, dims=2)[:]

    δdirtol = 1e-6 # too small directions get ignored ≈ 0 (solver precision issue) 
    sample_num = 0
    samplelength = 0
    updatesamplesizelength = true
    for n=1:N

        if updatesamplesizelength # switch to samples 
            direction_point = (@view wpoints[:, rand(1:nwpts)]) - (@view current_point[:]) # use warmup points to find direction in warmup phase
        else
            direction_point = (@view shat[:]) - (@view current_point[:]) # after warmup phase, only find directions in sampled space
        end

        λmax = 1e10
        λmin = -1e10 
        for i in eachindex(lbs)
            δlower = lbs[i] - current_point[i]
            δupper = ubs[i] - current_point[i]
            if direction_point[i] < -δdirtol
                lower = δupper/direction_point[i]
                upper = δlower/direction_point[i]
            elseif direction_point[i] > δdirtol
                lower = δlower/direction_point[i]
                upper = δupper/direction_point[i]
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
    
        λ = rand()*(λmax - λmin) + λmin # random step size
        current_point .= current_point .+ λ .* direction_point # will be feasible
        shat .= (n.*shat .+ current_point)./(n+1)

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

    return samples
end