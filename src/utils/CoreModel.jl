Base.isequal(model1::CoreModel, model2::CoreModel) =
    isequal(model1.S, model2.S) &&
    isequal(model1.b, model2.b) &&
    isequal(model1.c, model2.c) &&
    isequal(model1.xl, model2.xl) &&
    isequal(model1.xu, model2.xu) &&
    isequal(model1.rxns, model2.rxns) &&
    isequal(model1.mets, model2.mets)

Base.copy(model::CoreModel) =
    CoreModel(model.S, model.b, model.c, model.xl, model.xu, model.rxns, model.mets)

Base.isequal(model1::CoreCoupledModel, model2::CoreCoupledModel) =
    isequal(model1.lm, model2.lm) &&
    isequal(model1.C, model2.C) &&
    isequal(model1.cl, model2.cl) &&
    isequal(model1.cu, model2.cu)

Base.copy(model::CoreCoupledModel) = CoreCoupledModel(model.lm, model.C, model.cl, model.cu)
