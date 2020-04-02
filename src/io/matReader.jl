import ..LinearModel

using MAT


function loadModel(filePath::String, varName::String)

    file = matopen(filePath)
    vars = matread(filePath)

    if exists(file, varName)

        model = vars[varName]
        modelKeys = keys(model)

        S = model["S"]
        b = vec(model["b"])
        c = vec(model["c"])
        ub = vec(model["ub"])
        lb = vec(model["lb"])
        # rxns and mets get loaded as Array{Any,2} need to convert
        rxns = convert(Array{String, 1}, vec(model["rxns"]))
        mets = convert(Array{String, 1}, vec(model["mets"]))

        return LinearModel(S, b, c, lb, ub, rxns, mets)

    end
    close(file)
end

export loadModel
