Base.isequal(model1::LinearModel, model2::LinearModel) =
    model1.S == model2.S &&
    model1.b == model2.b &&
    model1.c == model2.c &&
    model1.xl == model2.xl &&
    model1.xu == model2.xu &&
    model1.rxns == model2.rxns &&
    model1.mets == model2.mets

Base.copy(model::LinearModel) =
    LinearModel(model.S, model.b, model.c, model.xl, model.xu, model.rxns, model.mets)

Base.isequal(model1::CoupledLinearModel, model2::CoupledLinearModel) =
    model1.lm == model2.lm &&
    model1.C == model2.C &&
    model1.cl == model2.cl &&
    model1.cu == model2.cu

Base.copy(model::CoupledLinearModel) =
    CoupledLinearModel(model.lm, model.C, model.cl, model.cu)
