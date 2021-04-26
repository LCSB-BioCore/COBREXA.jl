"""
    change_constraint(id::String, lb, ub)

Change the lower and upper bounds (`lb` and `ub` respectively) of reaction `id`.
"""
function change_constraint(id::String, lb, ub)
    (model, opt_model) -> begin
        ind = first(indexin([id], reactions(model)))
        set_bound(ind, opt_model, lb = lb, ub = ub)
    end
end

"""
    change_objective(objective_functions::Union{String,Vector{String}}; weights=[], sense=MOI.MAX_SENSE)

Callback function to change the objective function used in a constraint based analysis function. 
`objective_functions` can be a single reaction or an array of reactions (input their string `id`s).
Optionally specify their `weights` and the sense of the new objective (`MOI.MAX_SENSE`, `MOI.MIN_SENSE`).
Note, this function sets the sense of the objective to `MOI.MAX_SENSE` by default if not specified.
"""
function change_objective(
    objective_functions::Union{String,Vector{String}};
    weights = [],
    sense = MOI.MAX_SENSE,
)
    (model, opt_model) -> begin

        # Construct objective_indices array
        if typeof(objective_functions) == String
            objective_indices = [first(indexin([objective_functions], reactions(model)))]
        else
            objective_indices = [
                first(indexin([rxnid], reactions(model))) for rxnid in objective_functions
            ]
        end

        # Initialize weights
        opt_weights = spzeros(length(model.reactions))

        isempty(weights) && (weights = ones(length(objective_indices))) # equal weights

        for (j, i) in enumerate(objective_indices)
            opt_weights[i] = weights[j]
        end

        v = opt_model[:x]
        @objective(opt_model, sense, sum(opt_weights[i] * v[i] for i in objective_indices))
    end
end
