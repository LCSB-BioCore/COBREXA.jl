"""
$(TYPEDSIGNATURES)

Change the abundances of `model` with `new_abundances` in place.
"""
function change_abundances!(model::CommunityModel, new_abundances::Vector{Float64})
    isapprox(sum(new_abundances), 1.0; atol = constants.tolerance) ||
        throw(ArgumentError("The abundances do not sum to 1."))
    model.abundances .= new_abundances
    nothing
end

"""
$(TYPEDSIGNATURES)

Return a shallow copy of `model` with the abundances changed.
"""
function change_abundances(model::CommunityModel, new_abundances::Vector{Float64})
    m = copy(model)
    m.abundances = copy(model.abundances)
    m.abundances .= new_abundances
    m
end

"""
$(TYPEDSIGNATURES)

Change the environmental bound (environmental exchange) reaction `rid` in
`model`. If the associated bound is `nothing`, then it is not changed.
"""
change_environmental_bound!(
    model::CommunityModel,
    rid::String;
    lower_bound = nothing,
    upper_bound = nothing,
) = change_environmental_bounds!(
    model,
    [rid];
    lower_bounds = [lower_bound],
    upper_bounds = [upper_bound],
)

"""
$(TYPEDSIGNATURES)

Plural variant of [`change_environmental_bound!`](@ref).
"""
function change_environmental_bounds!(
    model::CommunityModel,
    rids::Vector{String};
    lower_bounds = fill(nothing, length(rids)),
    upper_bounds = fill(nothing, length(rids)),
)
    for (rid, lb, ub) in zip(rids, lower_bounds, upper_bounds)
        isnothing(lb) || model.environmental_exchange_reactions[rid][2] = lb
        isnothing(ub) || model.environmental_exchange_reactions[rid][3] = ub
    end
end

"""
$(TYPEDSIGNATURES)

Return a shallow copy of `model` with the environmental bound (environmental
exchange) reaction `rid` in `model` changed. If the associated bound is
`nothing`, then it is not changed.
"""
change_environmental_bound(
    model::CommunityModel,
    rid::String;
    lower_bound = nothing,
    upper_bound = nothing,
) = change_environmental_bounds(
    model,
    [rid];
    lower_bounds = [lower_bound],
    upper_bounds = [upper_bound],
)

"""
$(TYPEDSIGNATURES)

Plural variant of [`change_environmental_bound`](@ref).
"""
function change_environmental_bounds(
    model::CommunityModel,
    rids::Vector{String};
    lower_bounds = fill(nothing, length(rids)),
    upper_bounds = fill(nothing, length(rids)),
)
    m = copy(model)
    m.environmental_exchange_reactions = copy(model.environmental_exchange_reactions)
    for (k, (mid, _lb, _ub)) in model.environmental_exchange_reactions
        _idx = indexin([k], rids)
        if isnothing(_idx)
            m.environmental_exchange_reactions[k] = [mid, _lb, _ub]
        else
            idx = first(_idx)
            lb = isnothing(lower_bounds[idx]) ? _lb : lower_bounds[idx]
            ub = isnothing(upper_bounds[idx]) ? _ub : upper_bounds[idx]
            m.environmental_exchange_reactions[k] = [mid, lb, ub]
        end
    end
    m
end
