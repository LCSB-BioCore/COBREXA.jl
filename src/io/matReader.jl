import ..LinearModel

using MAT

"""
Load a model in MAT (Matlab) format and returns a `LinearModel`

See also: `MAT.jl`
"""
function loadModel(filePath::String, varName::String)

    # open file and read
    file = matopen(filePath)
    vars = matread(filePath)

    if exists(file, varName)

        model = vars[varName]
        modelKeys = ["S", "b", "c", "ub", "lb"]

        for key = modelKeys
            if !(key in keys(model))
                error("No variable $key found in the MAT file.")
            end
        end

        S = model["S"]
        b = vec(model["b"])
        c = vec(model["c"])
        ub = vec(model["ub"])
        lb = vec(model["lb"])
        # rxns and mets get loaded as Array{Any,2} need to convert
        rxns = conversionMetRxns(model["rxns"])
        mets = conversionMetRxns(model["mets"])

        return LinearModel(S, b, c, lb, ub, rxns, mets)

    else
        error("Variable `varName` does not exist in the specified MAT file.")
    end
    close(file)
end

function conversionMetRxns(x)
    if x isa Array{Any, 2}
        return convert(Array{String, 1}, vec(x))
    else
        return vec(x)
    end
end
