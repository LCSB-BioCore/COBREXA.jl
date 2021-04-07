"""
Convert a dictionary read from a MAT file to LinearModel
"""
function convert_to_linear_model(model::Dict)
    model_keys = ["S", "b", "c", "ub", "lb"]

    for key in model_keys
        if !(key in keys(model))
            error("No variable $key found in the MAT file.")
        end
    end

    S = model["S"]
    b = vec(model["b"])
    c = vec(model["c"])
    ub = vec(model["ub"])
    lb = vec(model["lb"])
    rxns = vec(string.(model["rxns"]))
    mets = vec(string.(model["mets"]))

    return LinearModel(S, b, c, lb, ub, rxns, mets)
end

"""
Convert a LinearModel to exportable format
SparseVectors are not written and read properly, SparseMatrix is okay
"""
function _convert_to_exportable(model::LinearModel)
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
