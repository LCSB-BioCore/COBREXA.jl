"""
Write a model into a MAT (Matlab) format

NB: Does NOT export general inequality constraints (eg coupling)

See also: `MAT.jl`
"""
function write_model(file_path::String, model::LinearModel, var_name::String = "model")
    matwrite(file_path, Dict(var_name => convert_to_exportable(model)))
end

"""
Convert a LinearModel to exportable format
SparseVectors are not written and read properly, SparseMatrix is okay
"""
function convert_to_exportable(model::LinearModel)
    return Dict(
        "S" => model.S,
        "b" => Array(model.b),
        "c" => Array(model.c),
        "ub" => Array(model.xu),
        "lb" => Array(model.xl),
        "rxns" => model.rxns,
        "mets" => model.mets,
    )
end
