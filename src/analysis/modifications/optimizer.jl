
"""
    change_sense(objective_sense)

Change the objective sense of optimization.
Possible arguments are `MOI.MAX_SENSE` and `MOI.MIN_SENSE`.

If you want to change the objective and sense at the same time, use
[`change_objective`](@ref) instead to do both at once.
"""
function change_sense(objective_sense)
    (model, opt_model) -> COBREXA.JuMP.set_objective_sense(opt_model, objective_sense)
end

"""
    change_optimizer(optimizer)

Change the JuMP optimizer used to run the optimization.

This may be used to try different approaches for reaching the optimum, and in
problems that may require different optimizers for different parts, such as the
[`parsimonious_flux_balance_analysis`](@ref).
"""
function change_optimizer(optimizer)
    (model, opt_model) -> COBREXA.JuMP.set_optimizer(opt_model, optimizer)
end

"""
    change_optimizer_attribute(attribute_key, value)

Change a JuMP optimizer attribute. The attributes are optimizer-specific, refer
to the JuMP documentation and the documentation of the specific optimizer for
usable keys and values.
"""
function change_optimizer_attribute(attribute_key, value)
    (model, opt_model) ->
        COBREXA.JuMP.set_optimizer_attribute(opt_model, attribute_key, value)
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
    change_objective(new_objective::Union{String,Vector{String}}; weights=[], sense=MOI.MAX_SENSE)

Modification that changes the objective function used in a constraint based
analysis function.  `new_objective` can be a single reaction identifier, or an
array of reactions identifiers.

Optionally, the objective can be weighted by a vector of `weights`, and a
optimization `sense` can be set.
"""
function change_objective(
    new_objective::Union{String,Vector{String}};
    weights = [],
    sense = MOI.MAX_SENSE,
)
    (model, opt_model) -> begin

        # Construct objective_indices array
        if typeof(new_objective) == String
            objective_indices = indexin([new_objective], reactions(model))
        else
            objective_indices =
                [first(indexin([rxnid], reactions(model))) for rxnid in new_objective]
        end

        any(isnothing.(objective_indices)) && throw(
            DomainError(new_objective, "No matching reaction found for one or more ids."),
        )

        # Initialize weights
        opt_weights = spzeros(n_reactions(model))

        isempty(weights) && (weights = ones(length(objective_indices))) # equal weights

        for (j, i) in enumerate(objective_indices)
            opt_weights[i] = weights[j]
        end

        v = opt_model[:x]
        @objective(opt_model, sense, sum(opt_weights[i] * v[i] for i in objective_indices))
    end
end
