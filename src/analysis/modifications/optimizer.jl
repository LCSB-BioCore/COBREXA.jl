
"""
    change_sense(objective_sense)

Change the objective sense of optimization.
Possible arguments are `MOI.MAX_SENSE` and `MOI.MIN_SENSE`.

If you want to change the objective and sense at the same time, use
[`change_objective`](@ref) instead to do both at once.
"""
function change_sense(objective_sense)
    (model, opt_model) -> set_objective_sense(opt_model, objective_sense)
end

"""
    change_optimizer(optimizer)

Change the JuMP optimizer used to run the optimization.

This may be used to try different approaches for reaching the optimum, and in
problems that may require different optimizers for different parts, such as the
[`parsimonious_flux_balance_analysis`](@ref).
"""
function change_optimizer(optimizer)
    (model, opt_model) -> set_optimizer(opt_model, optimizer)
end

"""
    change_optimizer_attribute(attribute_key, value)

Change a JuMP optimizer attribute. The attributes are optimizer-specific, refer
to the JuMP documentation and the documentation of the specific optimizer for
usable keys and values.
"""
function change_optimizer_attribute(attribute_key, value)
    (model, opt_model) -> set_optimizer_attribute(opt_model, attribute_key, value)
end

"""
    constrain_objective_value(tolerance)

Limit the objective value to `tolerance`-times the current objective value, as
with [`objective_bounds`](@ref).
"""
function constrain_objective_value(tolerance)
    return (model, opt_model) -> begin
        lambda_min, lambda_max = objective_bounds(tolerance)(objective_value(opt_model))
        old_objective = objective_function(opt_model)
        @constraint(opt_model, lambda_min <= sum(old_objective) <= lambda_max)
    end
end

"""
    change_constraint(id::String, lb, ub)

Change the lower and upper bounds (`lb` and `ub` respectively) of reaction `id`.
"""
function change_constraint(id::String, lb, ub)
    (model, opt_model) -> begin
        ind = first(indexin([id], reactions(model)))
        isnothing(ind) && throw(DomainError(id, "No matching reaction was found."))
        set_optmodel_bound!(ind, opt_model, lb = lb, ub = ub)
    end
end

"""
    change_objective(rxn::Union{Int,String}; weight::Float64 = 1.0, sense = MOI.MAX_SENSE)

Modification that changes the objective to maximize or minimize the specified
reaction (either by index or string ID), optionally specifying the objective
weight.
"""
change_objective(rxn::Union{Int,String}; weight::Float64 = 1.0, sense = MOI.MAX_SENSE) =
    change_objective([rxn]; weights = [weight], sense = sense)

"""
    change_objective(
        rids::Vector{String};
        weights::Vector{Float64} = ones(length(rids)),
        sense = MOI.MAX_SENSE,
    )

Modification that changes the objective function used in a constraint based
analysis function to maximize (or minimize, based on `sense`) the sum
(optionally weighted by `weights`) of the rates of the reactions specified by
string IDs.
"""
change_objective(
    rids::Vector{String};
    weights::Vector{Float64} = ones(length(rids)),
    sense = MOI.MAX_SENSE,
) =
    (model, opt_model) -> begin

        ridxs = indexin(rids, reactions(model))
        any(isnothing, ridxs) &&
            throw(DomainError(rids[isnothing.(ridxs)], "Unknown reaction IDs"))

        change_objective(Vector{Int}(ridxs); weights = weights, sense = sense)(
            model,
            opt_model,
        )
    end

"""
    change_objective(
        ridxs::Vector{Int};
        weights::Vector{Float64} = ones(length(ridxs)),
        sense = MOI.MAX_SENSE,
    )

A potentially more efficient variant of [`change_objective`](@ref) that works
on integer indexes of the reactions.
"""
change_objective(
    ridxs::Vector{Int};
    weights::Vector{Float64} = ones(length(ridxs)),
    sense = MOI.MAX_SENSE,
) =
    (model, opt_model) -> begin

        @objective(
            opt_model,
            sense,
            sum(weights[i] * opt_model[:x][ridx] for (i, ridx) in enumerate(ridxs))
        )
    end
