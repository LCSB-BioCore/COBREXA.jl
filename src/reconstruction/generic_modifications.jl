
"""
    change_sense(objective_sense)

Change the objective sense. 
Possible arguments are `MOI.MAX_SENSE` and `MOI.MIN_SENSE`.
Note, [`change_objective`](@ref) sets the sense of the objective,
so it doesn't make sense to use this function AND [`change_objective`](@ref) simultaneously. 
"""
function change_sense(objective_sense)
    (model, opt_model) -> begin
        COBREXA.JuMP.set_objective_sense(opt_model, objective_sense)
    end
end

"""
    change_solver(optimizer)

Change the solver (`optimizer`) used to solve the model.
Typically the solver is specified as a required argument in a function. 
However, this function is useful if the problem has multiple subparts that 
require different solvers.

See also: [`parsimonious_flux_balance_analysis`](@ref)
"""
function change_solver(optimizer)
    (model, opt_model) -> begin
        COBREXA.JuMP.set_optimizer(opt_model, optimizer)
    end
end

"""
    change_solver_attribute(option_key, option_val)

Change a solver attribute. These attributes are solver specific,
refer the either JuMP or the solver you are using's documentation.
"""
function change_solver_attribute(option_key, option_val)
    (model, opt_model) -> begin
        COBREXA.JuMP.set_optimizer_attribute(opt_model, option_key, option_val)
    end
end

"""
    constrain_objective_value(optimum_bound)

Limit the objective value to `optimum_bound`-times the current objective value.
"""
function constrain_objective_value(optimum_bound)
    return (model, opt_model) -> begin
        λ = objective_value(opt_model)
        λmin = min(optimum_bound * λ, λ * 1.0 / optimum_bound)
        λmax = max(optimum_bound * λ, λ * 1.0 / optimum_bound)
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
