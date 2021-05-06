"""
    _get_warmup_points(cbm; random_objective=false, numstop=1e10)

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
function _get_warmup_points(cbm; random_objective = false, numstop = 1e10)
    v = cbm[:x]
    ubs = cbm[:ubs]
    lbs = cbm[:lbs]
    # determine which rxns should be max/min-ized (non fixed fluxes)
    fixed_rxns = Int[]
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
