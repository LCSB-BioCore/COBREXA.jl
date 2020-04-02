import ..LinearModel

using MAT


function loadModel(filePath::String, varName::String)

    # open file and read
    file = matopen(filePath)
    vars = matread(filePath)

    if exists(file, varName)

        model = vars[varName]
        modelKeys = keys(model)

        if "S" in modelKeys
            S = model["S"]
        else
            error("No S matrix was found in the MAT file.")
        end

        if "b" in modelKeys
            b = vec(model["b"])
        else
            error("No b vector found in the MAT file.")
        end

        if "c" in modelKeys
            c = vec(model["c"])
        else
            error("No c vector found in the MAT file.")
        end

        if "ub" in modelKeys
            ub = vec(model["ub"])
        else
            error("No ub vector found in the MAT file.")
        end

        if "lb" in modelKeys
            lb = vec(model["lb"])
        else
            error("No lb vector found in the MAT file.")
        end

        # rxns and mets get loaded as Array{Any,2} need to convert

        if "rxns" in modelKeys
            rxns = convert(Array{String, 1}, vec(model["rxns"]))
        else
            error("No rxns vector found in the MAT file.")
        end

        if "mets" in modelKeys
            mets = convert(Array{String, 1}, vec(model["mets"]))
        else
            error("No mets vector found in the MAT file.")
        end

        return LinearModel(S, b, c, lb, ub, rxns, mets)

    else
        error("Variable `varName` does not exist in the specified MAT file.")
    end
    close(file)
end

export loadModel
