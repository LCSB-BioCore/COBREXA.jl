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
