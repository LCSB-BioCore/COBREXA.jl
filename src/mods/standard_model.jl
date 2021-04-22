"""
    change_constraint(reaction::Reaction, lb, ub)

Change the lower and upper bounds (`lb` and `ub` respectively) of `reaction`.
"""
function change_constraint(reaction::Reaction, lb, ub)
    (model, opt_model) -> begin
        ind = index_of(reaction.id, reactions(model))
        set_bound(ind, opt_model, lb = lb, ub = ub)
    end
end

"""
    change_objective(objective_functions::Union{Reaction, Vector{Reaction}}; weights=[], sense=MOI.MAX_SENSE)

Callback function to change the objective function used in a constraint based analysis function. 
`objective_functions` can be a single reaction or an array of reactions (of type `Reaction`).
Optionally specify their `weights` and the sense of the new objective (`MOI.MAX_SENSE`, `MOI.MIN_SENSE`).
Note, this function sets the sense of the objective to `MOI.MAX_SENSE` by default if not specified.
"""
function change_objective(
    objective_functions::Union{Reaction,Vector{Reaction}};
    weights = [],
    sense = MOI.MAX_SENSE,
)
    (model, opt_model) -> begin

        # Construct objective_indices array
        if typeof(objective_functions) == Reaction
            objective_indices = [index_of(objective_functions.id, reactions(model))]
        else
            objective_indices =
                [index_of(rxn.id, reactions(model)) for rxn in objective_functions]
        end

        # Initialize weights
        opt_weights = zeros(length(model.reactions))

        isempty(weights) && (weights = ones(length(objective_indices))) # equal weights

        wcounter = 1
        for i = 1:length(model.reactions)
            if i in objective_indices
                opt_weights[i] = weights[wcounter]
                wcounter += 1
            end
        end

        v = opt_model[:x]
        @objective(opt_model, sense, sum(opt_weights[i] * v[i] for i in objective_indices))
    end
end
