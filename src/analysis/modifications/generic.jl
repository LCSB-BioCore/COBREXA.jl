
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
    change_constraint(id::String; lb=nothing, ub=nothing)

Change the lower and upper bounds (`lb` and `ub` respectively) of reaction `id` if supplied.
"""
change_constraint(id::String; lb = nothing, ub = nothing) =
    (model, opt_model) -> begin
        ind = first(indexin([id], reactions(model)))
        isnothing(ind) && throw(DomainError(id, "No matching reaction was found."))
        set_optmodel_bound!(ind, opt_model, lb = lb, ub = ub)
    end

"""
    change_objective(new_objective::Union{String,Vector{String}}; weights=[], sense=MAX_SENSE)

Modification that changes the objective function used in a constraint based
analysis function.  `new_objective` can be a single reaction identifier, or an
array of reactions identifiers.

Optionally, the objective can be weighted by a vector of `weights`, and a
optimization `sense` can be set to either `MAX_SENSE` or `MIN_SENSE`.
"""
change_objective(
    new_objective::Union{String,Vector{String}};
    weights = [],
    sense = MAX_SENSE,
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

"""
    add_crowding_constraint(weights::Dict{Int64, Float64})

Adds a molecular crowding constraint to the optimization problem: `∑ wᵢ × vᵢ ≤ 1` where `wᵢ`
is a weight and `vᵢ` is a flux index in the model's reactions specified in `weights` as `vᵢ
=> wᵢ` pairs.

See Beg, Qasim K., et al. "Intracellular crowding defines the mode and sequence of
substrate uptake by Escherichia coli and constrains its metabolic activity." Proceedings of
the National Academy of Sciences 104.31 (2007) for more details.
"""
add_crowding_constraint(weights::Dict{Int64,Float64}) =
    (model, opt_model) -> begin
        idxs = collect(keys(weights)) # order of keys and values is the same
        ws = values(weights)
        # since fluxes can be positive or negative, need absolute value: ∑ wᵢ × |vᵢ| ≤ 1
        # introduce slack variables to handle this
        @variable(opt_model, crowding_slack[1:length(weights)])
        @constraint(opt_model, crowding_slack .>= opt_model[:x][idxs])
        @constraint(opt_model, crowding_slack .>= -opt_model[:x][idxs])
        @constraint(opt_model, sum(w * crowding_slack[i] for (i, w) in enumerate(ws)) <= 1)
    end

"""
    add_crowding_constraint(weight::Float64; kwargs)

Variant of [`add_crowding_constraint`](@ref) that takes a single weight and assigns it to
each internal reaction flux, where internal reactions are identified with
[`find_internal_reactions`](@ref) and `kwargs` are passed to this function.
"""
add_crowding_constraint(weight::Float64; kwargs...) =
    (model, opt_model) -> begin
        idxs = find_internal_reactions(model; kwargs...)
        add_crowding_constraint(Dict(zip(idxs, fill(weight, length(idxs)))))(
            model,
            opt_model,
        )
    end

"""
    add_crowding_constraint(weights::Dict{String, Float64})

Variant of [`add_crowding_constraint`](@ref) that takes a dictinary of reactions `ids`
instead of reaction indices mapped to weights.
"""
add_crowding_constraint(weights::Dict{String,Float64}) =
    (model, opt_model) -> begin
        idxs = indexin(keys(weights), reactions(model))
        nothing in idxs && throw(ArgumentError("Reaction id not found in model."))
        add_crowding_constraint(Dict(zip(Int.(idxs), values(weights))))(model, opt_model)
    end
