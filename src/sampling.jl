"""
get_warmup_points(cbmodel, v, mb, ubs, lbs)

Generate warmup points for all the reactions on the model that 
are not fixed. Assumes you feed in a JuMP model that is already
constrained by however you want it to be.
"""
function get_warmup_points(cbmodel, v, ubs, lbs; randomobjective=false)
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
    var_rxn_inds = filter(x->!(x in fixed_rxns), 1:length(v))

    wpoints = zeros(length(v), 2*length(var_rxn_inds))

    if randomobjective
        for i in 1:2:size(wpoints, 2)
            @objective(cbmodel, Max, sum(rand()*v[ii] for ii in var_rxn_inds))
            optimize!(cbmodel)
            for j=1:size(wpoints, 1)
                wpoints[j, i] = value(v[j])
            end
            
            @objective(cbmodel, Min, sum(rand()*v[ii] for ii in var_rxn_inds))
            optimize!(cbmodel)
            for j=1:size(wpoints, 1)
                wpoints[j, i+1] = value(v[j])
            end
        end
    else
        for (ii, i) in enumerate(1:2:size(wpoints, 2))
            @objective(cbmodel, Max, v[var_rxn_inds[ii]])
            optimize!(cbmodel)
            for j=1:size(wpoints, 1)
                wpoints[j, i] = value(v[j])
            end
            
            @objective(cbmodel, Min, v[var_rxn_inds[ii]])
            optimize!(cbmodel)
            for j=1:size(wpoints, 1)
                wpoints[j, i+1] = value(v[j])
            end
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

function achr(N::Int64, wpoints::Array{Float64, 2}, model::Model, ubcons, lbcons; burnin=0.9, keepevery=10)  
    S, _, _, _ = CobraTools.get_core_model(model) # assume S has not been modified from model
    ubs, lbs = CobraTools.get_bound_vectors(ubcons, lbcons)
    
    nwpts = size(wpoints, 2) # number of warmup points generated

    Nkeep = round(Int64, burnin*N) # start storing from here only 
    
    samples = zeros(size(wpoints, 1), round(Int64, length(Nkeep:N)/keepevery, RoundUp)) # sample storage

    center_point = mean(wpoints, dims=2)[:]

    w = zeros(size(wpoints, 1)) # direction vector
    λlower = zeros(size(S, 2))
    λupper = zeros(size(S, 2))

    δdirtol = 1e-6 # too small directions get ignored solver precision issue 

    current_point = zeros(size(wpoints, 1))
    current_point .= wpoints[:, rand(1:nwpts)] # initial point
    sample_num = 0

    for n=1:N
        foundit = false # found a feasible direction
        λmax = 0.0
        λmin = 0.0
        numiters = 0
        while numiters < 2*(nwpts+sample_num) # maximum time spent looking for feasible direction
            # pick a random direction from samples and warmup points
            if rand() < nwpts/(nwpts + sample_num)
                w .= wpoints[:, rand(1:nwpts)] .- center_point
            else
                w .= samples[:, rand(1:(sample_num))] .- center_point
            end
            w .= w./norm(w) # direction

            for i in eachindex(w)
                δlower = lbs[i] - current_point[i]
                δupper = ubs[i] - current_point[i]
                if w[i] < -δdirtol
                    λlower[i] = δupper/w[i]
                    λupper[i] = δlower/w[i]    
                elseif w[i] > δdirtol
                    λlower[i] = δlower/w[i]
                    λupper[i] = δupper/w[i]
                else
                    λlower[i] = -1e3
                    λupper[i] = 1e3
                end
            end

            λmax = minimum(λupper)
            λmin = maximum(λlower)
            
            if λmax > λmin
                foundit = true
                break
            else
                numiters += 1
            end
        end

        if foundit == false
            @warn "Error:  no feasible direction found."
            break
        end

        λ = rand()*(λmax - λmin) + λmin
        current_point .= current_point .+ λ .* w

        center_point .= ((nwpts + n - 1).*center_point .+ current_point) ./ (nwpts + n)

        if n >= Nkeep && n % keepevery == 0
            sample_num += 1
            samples[:, sample_num] .= current_point
        end
        
    end

    return samples
end

