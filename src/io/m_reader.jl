"""
    read_matlab_model(file_location::String)
"""
function read_matlab_model(file_location::String)
    matfile = matread(file_location)
    model_name = collect(keys(matfile))[1]
    modeldict = matfile[model_name]

    # the model_id can be written in many places, try varying levels of specificity
    model_id = haskey(modeldict, "description") ? modeldict["description"] : model_name
    model_id = haskey(modeldict, "modelName") ? modeldict["modelName"] : model_name # more specific

    mets = Metabolite[]
    for i in eachindex(modeldict["mets"])
        id = haskey(modeldict, "mets") ? modeldict["mets"][i] : ""
        if id == ""
            continue
        end

        name = haskey(modeldict, "metNames") ? modeldict["metNames"][i] : ""
        compartment = ""
        formula = ""
        if haskey(modeldict, "metFormulas")
            formula = string(modeldict["metFormulas"][i])
        elseif haskey(modeldict, "metFormula")
            formula = string(modeldict["metFormula"][i])
        end

        charge = 0 # sometimes inconsistently named
        if haskey(modeldict, "metCharge") && !isnan(modeldict["metCharge"][i])
            charge = modeldict["metCharge"][i]
        elseif haskey(modeldict, "metCharges") && !isnan(modeldict["metCharges"][i])
            charge = modeldict["metCharges"][i]
        end

        # look for annotation data, assume delimited by "; "
        annotation = Dict{String,Union{Array{String,1},String}}()
        anno_kid = Dict(
            "metBiGGID" => "bigg.metabolite",
            "metKEGGID" => "kegg.compound",
            "metMetaNetXID" => "metanetx.chemical",
            "metChEBIID" => "chebi",
        )
        for (anno, kid) in anno_kid
            if haskey(modeldict, anno)
                annotation[kid] = string.(split(string(modeldict[anno][i]), "; "))
            end
        end
        if haskey(modeldict, "metSBOTerms")
            annotation["sbo"] = string(modeldict["metSBOTerms"][i])
        end

        # look for some notes
        notes = Dict{String,Array{String,1}}()
        if haskey(modeldict, "metNotes")
            notes["note"] = string.(split(string(modeldict["metNotes"][i]), "; "))
        end

        push!(mets, Metabolite(id, name, formula, charge, compartment, notes, annotation))
    end

    genes = Gene[]
    for i in eachindex(modeldict["genes"])
        id = haskey(modeldict, "genes") ? modeldict["genes"][i] : ""
        if id == ""
            continue # skip blanks
        end

        # these fields often don't exist in the matlab models, ignore for now
        name = ""
        notes = Dict{String,Array{String,1}}()
        annotation = Dict{String,Union{Array{String,1},String}}()

        push!(genes, Gene(id, name, notes, annotation))
    end

    rxns = Reaction[]
    for i in eachindex(modeldict["rxns"])
        id = haskey(modeldict, "rxns") ? modeldict["rxns"][i] : ""
        if id == ""
            continue # skip blanks
        end

        name = haskey(modeldict, "rxnNames") ? modeldict["rxnNames"][i] : ""
        metinds = findall(x -> x .!= 0.0, modeldict["S"][:, i])
        metabolites =
            Dict{Metabolite,Float64}(mets[j] => modeldict["S"][j, i] for j in metinds)

        lb = haskey(modeldict, "lb") ? modeldict["lb"][i] : -1000.0 # reversible by default
        ub = haskey(modeldict, "ub") ? modeldict["ub"][i] : 1000.0 # reversible by default

        grr_string = haskey(modeldict, "grRules") ? modeldict["grRules"][i] : ""
        subsystem = join(modeldict["subSystems"][i], "; ")

        objective_coefficient = haskey(modeldict, "c") ? modeldict["c"][i] : 0.0

        # look for some annotations
        annotation = Dict{String,Union{Array{String,1},String}}()
        anno_kids = Dict(
            "rxnKEGGID" => "kegg.reaction",
            "rxnECNumbers" => "ec-code",
            "rxnBiGGID" => "bigg.reaction",
        )
        for (anno, kid) in anno_kids
            if haskey(modeldict, anno)
                annotation[kid] = string.(split(string(modeldict[anno][i]), "; "))
            end
        end
        if haskey(modeldict, "rxnSBOTerms")
            annotation["sbo"] = string(modeldict["rxnSBOTerms"][i])
        end

        # look for some notes
        notes = Dict{String,Array{String,1}}()
        if haskey(modeldict, "rxnNotes")
            notes["note"] = string.(split(string(modeldict["rxnNotes"][i]), "; "))
        end

        # get gene reaction rule
        grr = parse_grr(grr_string, genes)

        push!(
            rxns,
            Reaction(
                id,
                name,
                metabolites,
                lb,
                ub,
                grr,
                subsystem,
                notes,
                annotation,
                objective_coefficient,
            ),
        )
    end

    return CobraModel(model_id, rxns, mets, genes)
end
