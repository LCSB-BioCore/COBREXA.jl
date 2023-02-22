"""
$(TYPEDSIGNATURES)

Change the abundances of `model` with `new_abundances` in place.
"""
function change_abundances!(model::CommunityModel, new_abundances::Vector{Float64})
    check_abundances(new_abundances)
    model.abundances .= new_abundances
    nothing
end

"""
$(TYPEDSIGNATURES)

Return a shallow copy of `model` with the abundances changed.
"""
function change_abundances(model::CommunityModel, new_abundances::Vector{Float64})
    check_abundances(new_abundances)
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
    idxs = check_environmental_ids(model, rids)
    for (idx, lb, ub) in zip(idxs, lower_bounds, upper_bounds)
        isnothing(lb) || (model.environmental_links[idx].lower_bound = lb)
        isnothing(ub) || (model.environmental_links[idx].upper_bound = ub)
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
    idxs = check_environmental_ids(model, rids)
    m = copy(model)
    m.environmental_links = copy(model.environmental_links)
    for (idx, lb, ub) in zip(idxs, lower_bounds, upper_bounds)
        m.environmental_links[idx] = copy(model.environmental_links[idx])
        m.environmental_links[idx].reaction_id = model.environmental_links[idx].reaction_id
        m.environmental_links[idx].metabolite_id =
            model.environmental_links[idx].metabolite_id
        m.environmental_links[idx].lower_bound =
            isnothing(lb) ? model.environmental_links[idx].lower_bound : lb
        m.environmental_links[idx].upper_bound =
            isnothing(ub) ? model.environmental_links[idx].upper_bound : ub
    end
    m
end

"""
$(TYPEDSIGNATURES)

Return an [`EqualGrowthCommunityModel`](@ref) wrapper around `model`, optionally
specifying the `community_objective_id`.
"""
make_EqualGrowthCommunityModel(
    model::CommunityModel;
    community_objective_id = "equal_growth_rates_biomass_function",
) = EqualGrowthCommunityModel(model, community_objective_id)
