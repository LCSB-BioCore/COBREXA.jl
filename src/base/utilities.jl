Base.isequal(model1::LinearModel, model2::LinearModel) =
    isequal(model1.S, model2.S) &&
    isequal(model1.b, model2.b) &&
    isequal(model1.c, model2.c) &&
    isequal(model1.xl, model2.xl) &&
    isequal(model1.xu, model2.xu) &&
    isequal(model1.rxns, model2.rxns) &&
    isequal(model1.mets, model2.mets)

Base.copy(model::LinearModel) =
    LinearModel(model.S, model.b, model.c, model.xl, model.xu, model.rxns, model.mets)

Base.isequal(model1::CoupledLinearModel, model2::CoupledLinearModel) =
    isequal(model1.lm, model2.lm) &&
    isequal(model1.C, model2.C) &&
    isequal(model1.cl, model2.cl) &&
    isequal(model1.cu, model2.cu)

Base.copy(model::CoupledLinearModel) =
    CoupledLinearModel(model.lm, model.C, model.cl, model.cu)


"""
    set_bound(index, optimization_model; ub=1000, lb=-1000)
Helper function to set the bounds of variables.
The JuMP `set_normalized_rhs` function is a little confusing, so this function simplifies setting
constraints. Just supply the constraint `index` and the model and they will be set to `ub` and `lb`.
"""
function set_bound(vind, opt_model; ub = 1000, lb = -1000)
    if lb <= 0
        set_normalized_rhs(opt_model[:lbs][vind], abs(lb))
    else
        set_normalized_rhs(opt_model[:lbs][vind], -abs(lb))
    end
    set_normalized_rhs(opt_model[:ubs][vind], ub)
end

"""
    modify_constraint(reaction::Reaction, lb, ub)

Modify constraints of model reaction.
"""
function modify_constraint(reaction::Reaction, lb, ub)
    (model, opt_model) -> begin
        ind = model.reactions[reaction]
        set_bound(ind, opt_model, lb=lb, ub=ub)
    end
end

"""
    modify_sense(sense)

    Modify the objective sense. 
Allowed options are MOI.MAX_SENSE and MOI.MIN_SENSE.
"""
function modify_sense(sense)
    (model, opt_model) -> begin
        set_objective_sense(opt_model, sense)
    end
end

"""
modify_solver_attributes(option_key, option_value)

"""
function modify_solver_attributes(option_key, option_value)
    (model, opt_model) -> begin
        set_optimizer_attribute(opt_model, option_key, option_value)
    end
end

"""
    modify_objective(objective_function::Union{Reaction,Array{Reaction,1}}; objective_weights=[], sense=MOI.MAX_SENSE)

Modify the objective function and optionally specify weights 
(useful if more than one reaction is set in `objective_rxns`)
and the sense of the objective.

# Example
"""
function modify_objective(objective_function::Union{Reaction,Array{Reaction,1}}; objective_weights=[], sense=MOI.MAX_SENSE)
    (model, opt_model) -> begin
        v = opt_model[:x]

        # ensure that an array of objective indices are fed in
        if typeof(objective_function) == Reaction
            objective_indices = [model[objective_function]]
        else
            objective_indices = [model[rxn] for rxn in objective_function]
        end

        if isempty(objective_weights)
            weights = ones(length(objective_indices))
        end
        opt_weights = zeros(length(model.reactions))

        # update the objective function tracker
        wcounter = 1
        for i in eachindex(model.reactions)
            if i in objective_indices
                opt_weights[i] = weights[wcounter]
                wcounter += 1
            end
        end
        
        @objective(cbm, MOI.MAX_SENSE, sum(opt_weights[i] * v[i] for i in objective_indices))
    end
end
