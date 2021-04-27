
"""
    load_mat_model(file_name::String)

Load and return a MATLAB file `file_name` that contains a COBRA-compatible
model.
"""
function load_mat_model(file_name::String)::MATModel
    model_pair = first(matread(file_name))
    @_io_log @info "Loading MAT: taking a model with ID $(model_pair.first)"
    return MATModel(model_pair.second)
end

#TODO: this needs to get merged into MAT->StdModel conversion
function Base.convert(::Type{StandardModel}, mm::MATModel)
    modeldict = mm.mat

    # the model_id can be written in many places, try varying levels of specificity
    model_id = haskey(modeldict, "description") ? modeldict["description"] : model_name
    model_id = haskey(modeldict, "modelName") ? modeldict["modelName"] : model_name # more specific

    mets = OrderedDict{String,Metabolite}()
    for i in eachindex(modeldict["mets"])
        met = Metabolite()

        if haskey(modeldict, "mets")
            met.id = modeldict["mets"][i]
        else
            continue
        end

        if haskey(modeldict, "metNames")
            met.name = modeldict["metNames"][i]
        end

        if haskey(modeldict, "metFormulas")
            met.formula = string(modeldict["metFormulas"][i])
        end

        if haskey(modeldict, "metFormula")
            met.formula = string(modeldict["metFormula"][i])
        end

        if haskey(modeldict, "metCharge") && !isnan(modeldict["metCharge"][i])
            met.charge = modeldict["metCharge"][i]
        end

        if haskey(modeldict, "metCharges") && !isnan(modeldict["metCharges"][i])
            met.charge = modeldict["metCharges"][i]
        end

        # look for annotation data, assume delimited by "; "
        anno_kid = Dict(
            "metBiGGID" => "bigg.metabolite",
            "metKEGGID" => "kegg.compound",
            "metMetaNetXID" => "metanetx.chemical",
            "metChEBIID" => "chebi",
        )
        for (anno, kid) in anno_kid
            if haskey(modeldict, anno)
                met.annotation[kid] = string.(split(string(modeldict[anno][i]), "; "))
            end
        end

        if haskey(modeldict, "metSBOTerms")
            met.annotation["sbo"] = [string(modeldict["metSBOTerms"][i])]
        end

        if haskey(modeldict, "metNotes")
            met.notes["note"] = string.(split(string(modeldict["metNotes"][i]), "; "))
        end

        mets[met.id] = met
    end

    genes = OrderedDict{String,Gene}()
    for i in eachindex(modeldict["genes"])
        gene = Gene()
        if haskey(modeldict, "genes")
            gene.id = modeldict["genes"][i]
        else
            continue
        end

        # name and other fields don't generally exist in the matlab models,
        # ignoring them for now

        genes[gene.id] = gene
    end

    rxns = OrderedDict{String,Reaction}()
    for i in eachindex(modeldict["rxns"])
        rxn = Reaction()
        rxn.id = modeldict["rxns"][i]

        if haskey(modeldict, "rxnNames")
            rxn.name = modeldict["rxnNames"][i]
        end

        metinds = findall(x -> x .!= 0.0, modeldict["S"][:, i])
        rxn.metabolites = Dict{String,Float64}(
            mets[sj].id => modeldict["S"][j, i] for
            (sj, j) in zip(modeldict["mets"][metinds], metinds)
        )

        if haskey(modeldict, "lb")
            rxn.lb = modeldict["lb"][i]
        end

        if haskey(modeldict, "ub")
            rxn.ub = modeldict["ub"][i]
        end

        if haskey(modeldict, "grRules")
            grr_string = modeldict["grRules"][i]
            rxn.grr = _parse_grr(grr_string)
        end

        rxn.subsystem = join(modeldict["subSystems"][i], "; ")

        if haskey(modeldict, "c")
            rxn.objective_coefficient = modeldict["c"][i]
        end

        # look for some annotations
        anno_kids = Dict(
            "rxnKEGGID" => "kegg.reaction",
            "rxnECNumbers" => "ec-code",
            "rxnBiGGID" => "bigg.reaction",
        )

        for (anno, kid) in anno_kids
            if haskey(modeldict, anno)
                rxn.annotation[kid] = string.(split(string(modeldict[anno][i]), "; "))
            end
        end

        if haskey(modeldict, "rxnSBOTerms")
            rxn.annotation["sbo"] = [string(modeldict["rxnSBOTerms"][i])]
        end

        # look for some notes
        if haskey(modeldict, "rxnNotes")
            rxn.notes["note"] = string.(split(string(modeldict["rxnNotes"][i]), "; "))
        end

        rxns[rxn.id] = rxn
    end

    return StandardModel(model_id; reactions = rxns, metabolites = mets, genes = genes)
end
