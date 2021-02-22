"""
Write a model into a MAT (Matlab) format

NB: Does NOT export general inequality constraints (eg coupling)

See also: `MAT.jl`
"""
function writeModel(filePath::String, model::LinearModel, varName::String="model")
    matwrite(filePath, Dict(varName => convertToExportable(model)))
end

"""
Convert a LinearModel to exportable format
SparseVectors are not written and read properly, SparseMatrix is okay
"""
function convertToExportable(model::LinearModel)
    return Dict("S"    => model.S,
                "b"    => Array(model.b),
                "c"    => Array(model.c),
                "ub"   => Array(model.xu),
                "lb"   => Array(model.xl),
                "rxns" => model.rxns,
                "mets" => model.mets)
end
