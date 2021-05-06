
"""
    change_sense(objective_sense)

Change the objective sense of optimization.
Possible arguments are `MOI.MAX_SENSE` and `MOI.MIN_SENSE`.

If you want to change the objective and sense at the same time, use
[`change_objective`](@ref) instead to do both at once.
"""
function change_sense(objective_sense)
    (model, opt_model) ->
        COBREXA.JuMP.set_objective_sense(opt_model, objective_sense)
end

"""
    change_solver(optimizer)

Change the JuMP solver used to run the optimization.

This may be used to try different approaches for reaching the optimum, and in
problems that may require different solvers for different parts, such as the
[`parsimonious_flux_balance_analysis`](@ref).
"""
function change_solver(optimizer)
    (model, opt_model) ->
        COBREXA.JuMP.set_optimizer(opt_model, optimizer)
end

"""
    change_solver_attribute(option_key, option_val)

Change a JuMP solver attribute. These attributes are solver specific,
refer the either JuMP or the solver you are using's documentation.
"""
function change_solver_attribute(option_key, option_val)
    (model, opt_model) ->
        COBREXA.JuMP.set_optimizer_attribute(opt_model, option_key, option_val)
end

"""
    constrain_objective_value(tolerance)

Limit the objective value to `tolerance`-times the current objective value, as
with [`objective_bounds`](@ref).
"""
function constrain_objective_value(tolerance)
    return (model, opt_model) -> begin
        λmin, λmax = objective_bounds(tolerance)(objective_value(opt_model))
        old_objective = objective_function(opt_model)
        @constraint(opt_model, λmin <= sum(old_objective) <= λmax)
    end
end

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
            objective_indices = [first(indexin([new_objective], reactions(model)))]
        else
            objective_indices = [
                first(indexin([rxnid], reactions(model))) for rxnid in new_objective
            ]
        end

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


"""
    knockout(gene_ids::Array{String,1})

Callback function to set bounds of all reactions to zero which are affected by knocking out their respective genes
"""
function knockout(gene_ids::Vector{String})
    # Dev note: the three nested for loops are inefficiency. However:
    # - gene_ids (user input) will be probably only very few items
    # - model.genes[gene_id].reactions are just a few reactions (most genes don't code for a lot of reactions)
    # - reaction.grr also should only hold few items (reactions aren't coded by many different combinations of genes)
    # Let's avoid premature optimization for now and see if anyone ever has problems with this
    return (model, opt_model) -> begin
        all_reactions = reactions(model)
        for gene_id in gene_ids
            for reaction_id in gene_associated_reactions(model, gene_id)
                if all(
                    any(occursin.(gene_ids, gene_array)) > 0 for
                    gene_array in reaction_gene_association(model, reaction_id)
                )
                    set_bound(
                        first(indexin([reaction_id], all_reactions)),
                        opt_model,
                        ub = 0,
                        lb = 0,
                    )
                    # Also set bounds for model object
                    set_bound(model, reaction_id, ub = 0, lb = 0)
                end
            end
        end
    end
end

"""
    knockout(gene_id::String)

Callback function to set bounds of a reaction to zero which are affected by knocking out their respective genes
"""
function knockout(gene_id::String)
    return knockout([gene_id])
end
