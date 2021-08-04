
"""
    constrain_objective_value(tolerance)

Limit the objective value to `tolerance`-times the current objective value, as
with [`objective_bounds`](@ref).
"""
constrain_objective_value(tolerance) =
    (_, opt_model) -> begin
        lambda_min, lambda_max = objective_bounds(tolerance)(objective_value(opt_model))
        old_objective = objective_function(opt_model)
        @constraint(opt_model, lambda_min <= sum(old_objective) <= lambda_max)
    end

"""
    change_constraint(id::String, lb, ub)

Change the lower and upper bounds (`lb` and `ub` respectively) of reaction `id`.
"""
change_constraint(id::String, lb, ub) =
    (model, opt_model) -> begin
        ind = first(indexin([id], reactions(model)))
        isnothing(ind) && throw(DomainError(id, "No matching reaction was found."))
        set_optmodel_bound!(ind, opt_model, lb = lb, ub = ub)
    end

"""
    change_objective(new_objective::Union{String,Vector{String}}; weights=[], sense=MOI.MAX_SENSE)

Modification that changes the objective function used in a constraint based
analysis function.  `new_objective` can be a single reaction identifier, or an
array of reactions identifiers.

Optionally, the objective can be weighted by a vector of `weights`, and a
optimization `sense` can be set.
"""
change_objective(
    new_objective::Union{String,Vector{String}};
    weights = [],
    sense = MOI.MAX_SENSE,
) =
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
