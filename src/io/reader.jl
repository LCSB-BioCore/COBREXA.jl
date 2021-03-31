"""
Load a model in MAT (Matlab) format and returns a `LinearModel`

See also: `MAT.jl`
"""
function load_model(file_path::String, var_name::String)

    # read file
    vars = matread(file_path)

    if haskey(vars, var_name)
        return convert_to_linear_model(vars[var_name])
    else
        error("Variable `$var_name` does not exist in the specified MAT file.")
    end
end

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
