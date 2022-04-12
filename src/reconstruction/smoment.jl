"""
    change_bound!(model::SMomentModel, id; lower=nothing, upper=nothing)

Change the bound of variable in `model`. Does not change the bound if respective
bound is `nothing`. Note, for `SMomentModel`s, if the model used to construct the
`SMomentModel` has irreversible reactions, then these reactions will be
permanently irreversible in the model, i.e. changing their bounds to make them
reversible will have no effect.
"""
function change_bound!(model::SMomentModel, id; lower = nothing, upper = nothing)


    flux_for_idx =
        findfirst(x -> x == id * "§ARM§FOR" || x == id * "§FOR", model.irrev_reaction_ids)
    if !isnothing(flux_for_idx)
        if !isnothing(lower)
            if lower <= 0
                model.xl[flux_for_idx] = 0
            else
                model.xl[flux_for_idx] = lower
            end
        end
        if !isnothing(upper)
            if upper <= 0
                model.xu[flux_for_idx] = 0
            else
                model.xu[flux_for_idx] = upper
            end
        end
    end

    flux_rev_idx =
        findfirst(x -> x == id * "§ARM§REV" || x == id * "§REV", model.irrev_reaction_ids)
    if !isnothing(flux_rev_idx)
        if !isnothing(lower)
            if lower >= 0
                model.xu[flux_rev_idx] = 0
            else
                model.xu[flux_rev_idx] = -lower
            end
            if !isnothing(upper)
                if upper >= 0
                    model.xl[flux_rev_idx] = 0
                else
                    model.xl[flux_rev_idx] = -upper
                end
            end
        end
    end

    return nothing
end

"""
    change_bounds!(model::SMomentModel, ids; lower=fill(nothing, length(ids)), upper=fill(nothing, length(ids)))

Change the bounds of multiple variables in `model` simultaneously. See
[`change_bound`](@ref) for details.
"""
function change_bounds!(
    model::SMomentModel,
    ids;
    lower = fill(nothing, length(ids)),
    upper = fill(nothing, length(ids)),
)
    for (id, lower, upper) in zip(ids, lower, upper)
        change_bound!(model, id; lower = lower, upper = upper)
    end
end
