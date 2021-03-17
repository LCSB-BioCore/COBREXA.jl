"""
`n` by `n` sparse identity matrix.
"""
speye(n) = sparse(1.0I, n, n)

function Base.isequal(model1::LinearModel, model2::LinearModel)
    return model1.S == model2.S &&
           model1.b == model2.b &&
           model1.C == model2.C &&
           model1.cl == model2.cl &&
           model1.cu == model2.cu &&
           model1.c == model2.c &&
           model1.xl == model2.xl &&
           model1.xu == model2.xu &&
           model1.rxns == model2.rxns &&
           model1.mets == model2.mets
end

function Base.copy(model::LinearModel)
    return LinearModel(
        model.S,
        model.b,
        model.C,
        model.cl,
        model.cu,
        model.c,
        model.xl,
        model.xu,
        model.rxns,
        model.mets,
    )
end
