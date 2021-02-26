"""
    read_model(file_location::String))

Reads a model at file_location and returns a constraint based model (::CobraTools.Model).
Currently supported formats include SBML (.xml), Matlab (.mat) and JSON (.json) models.
The model format is inferred from the file extension.

Note, some meta-information may be lost when importing a model. Importantly, only information regarding the
reactions, metabolites and genes are imported. Currently reading JSON models captures the most meta-information
regarding reactions, metabolites and genes (e.g. the notes and annotation fields). 

When importing Matlab models some annotation and notes may not be imported because of the non-standard field names used by some models.
Gene reaction rules are successfully imported only if they adhere to this format: `"(YIL010W and YLR043C) or (YIL010W and YGR209C)"`. Other gene reaction rules
formats are not supported yet. 

In all cases the basic information should be imported, e.g. stoichiometrix matrix, constraints etc..
Advanced tools that require, e.g. metabolite formulas, gene reaction rules, and KEGG or BIGG IDs, will not function if these are improperly imported.
Always inspect the imported model before running analysis. 
"""
function read_model(file_location::String)
    if endswith(file_location, ".json")
        try 
            model = reconstruct_model_json(file_location)
        catch err
            @error "JSON model reading error.\n$err"
            model = CobraTools.Model()
        end
    elseif endswith(file_location, ".xml")
        try
            model = reconstruct_model_sbml(file_location)
        catch err
            @error "SBML model reading error.\n$err"
            model = CobraTools.Model()
        end
    elseif endswith(file_location, ".mat")
       try
            model = reconstruct_model_matlab(file_location)
       catch err
            @error "Matlab model reading error.\n$err"
            model = CobraTools.Model()
       end
    else
        @error "Model format not supported. The format is inferred from the file extension. Supported formats: *.mat, *.xml, *.json."
        model = CobraTools.Model()
    end
    return model
end

"""
    reconstructmodeljson(modeldict::String)
"""
function reconstruct_model_json(file_location::String)
    modeldict = JSON.parsefile(file_location)
    id = modeldict["id"]

    mets = Metabolite[]
    for met in modeldict["metabolites"]
        m = Metabolite(met)
        push!(mets, m)
    end
    
    rxns = Reaction[]
    for rxn in modeldict["reactions"]
        r = Reaction(rxn, mets)
        push!(rxns, r)
    end

    genes = Gene[]
    for gene in modeldict["genes"]
        g = Gene(gene)
        push!(genes, g)
    end

    grrs = Dict{String,  Array{Array{String, 1}, 1}}()
    for rxn in rxns
        isempty(rxn.grr) ? continue : (grrs[rxn.id] = parse_grr(rxn.grr))
    end
    
    return CobraTools.Model(id, rxns, mets, genes, grrs)
end

"""
    reconstruct_model_matlab(file_location::String)
"""
function reconstruct_model_matlab(file_location::String)
    matfile = matread(file_location)
    model_name = collect(keys(matfile))[1]
    modeldict = matfile[model_name]

    model_id = haskey(modeldict, "description") ? modeldict["description"] : model_name
    model_id = haskey(modeldict, "modelName") ? modeldict["modelName"] : model_name # more specific
    
    mets = Metabolite[]
    for i in eachindex(modeldict["mets"])
        id = haskey(modeldict, "mets") ? modeldict["mets"][i] : ""
        name = haskey(modeldict, "metNames") ? modeldict["metNames"][i] : ""
        compartment = ""
        formula = ""
        if haskey(modeldict, "metFormulas") 
            formula = modeldict["metFormulas"][i]
        elseif haskey(modeldict, "metFormula") 
            formula = modeldict["metFormula"][i]    
        end

        charge = 0 # sometimes inconsistently named
        if haskey(modeldict, "metCharge") && !isnan(modeldict["metCharge"][i])
            charge = modeldict["metCharge"][i]
        elseif haskey(modeldict, "metCharges") && !isnan(modeldict["metCharges"][i])
            charge = modeldict["metCharges"][i]
        end
        
        annotation = Dict{String, Union{Array{String, 1}, String}}()
        if haskey(modeldict, "metBiGGID")
            annotation["bigg.metabolite"] = [modeldict["metBiGGID"][i]]
        end
        if haskey(modeldict, "metSBOTerms")
            annotation["sbo"] = modeldict["metSBOTerms"][i]
        end
        if haskey(modeldict, "metKEGGID")
            annotation["kegg.compound"] = [modeldict["metKEGGID"][i]]
        end
        if haskey(modeldict, "metMetaNetXID")
            annotation["metanetx.chemical"] = [modeldict["metMetaNetXID"][i]]
        end
        if haskey(modeldict, "metChEBIID")
            annotation["chebi"] = [modeldict["metChEBIID"][i]]
        end
        
        notes = Dict{String, Array{String, 1}}()
        if haskey(modeldict, "metNotes")
            notes["note"] = [modeldict["metNotes"][i]]
        end            

        push!(mets, Metabolite(id, name, formula, charge, compartment, notes, annotation))
    end

    rxns = Reaction[]
    for i in eachindex(modeldict["rxns"])
        id = haskey(modeldict, "rxns") ? modeldict["rxns"][i] : ""
        name = haskey(modeldict, "rxnNames") ? modeldict["rxnNames"][i] : ""
        metinds = findall(x -> x .!= 0.0, modeldict["S"][:, i])
        metabolites = Dict{Metabolite, Float64}(mets[j]=>modeldict["S"][j, i] for j in metinds)
        lb = haskey(modeldict, "lb") ? modeldict["lb"][i] : -1000.0 # reversible by default
        ub = haskey(modeldict, "ub") ? modeldict["ub"][i] : 1000.0 # reversible by default
        grr = haskey(modeldict, "grRules") ? modeldict["grRules"][i] : ""
        subsystem = join(modeldict["subSystems"][i], "; ")
        objective_coefficient = 0.0
        objective_coefficient = haskey(modeldict, "c") ? modeldict["c"][i] : 0.0

        annotation = Dict{String, Union{Array{String, 1}, String}}()
        if haskey(modeldict, "rxnKEGGID")
            annotation["kegg.reaction"] = [modeldict["rxnKEGGID"][i]]
        end
        if haskey(modeldict, "rxnECNumbers")
            annotation["ec-code"] = string.(split(modeldict["rxnECNumbers"][i], "; "))
        end
        if haskey(modeldict, "rxnBiGGID")
            annotation["bigg.reaction"] = [modeldict["rxnBiGGID"][i]]
        end     
        
        notes = Dict{String, Array{String, 1}}()
        if haskey(modeldict, "rxnNotes")
            notes["note"] = [modeldict["rxnNotes"][i]]
        end

        push!(rxns, Reaction(id, name, metabolites, lb, ub, grr, subsystem, notes, annotation, objective_coefficient))
    end

    genes = Gene[]
    for i in eachindex(modeldict["genes"])
        id = haskey(modeldict, "genes") ? modeldict["genes"][i] : ""
        
        # these fields often don't exist in the matlab models
        name = ""
        notes =  Dict{String, Array{String, 1}}()
        annotation = Dict{String, Union{Array{String, 1}, String}}()
        
        push!(genes, Gene(id, name, notes, annotation))
    end

    grrs = Dict{String,  Array{Array{String, 1}, 1}}()
    for (i, rxn) in enumerate(rxns)
        if haskey(modeldict, "grRules") && !isempty(modeldict["grRules"][i])
            grrs[rxn.id] = parse_grr(modeldict["grRules"][i])
        end
    end
    
    return CobraTools.Model(model_id, rxns, mets, genes, grrs)
end

"""
    reconstruct_model_sbml(file_location::String)
"""
function reconstruct_model_sbml(file_location::String)
    m = readSBML(file_location)

# m is now a Model structure with:
# m.reactions
# m.species
# m.compartments
    # return Model()
    return m
end

"""
    parsegrr(string_rule)

Format: (YIL010W and YLR043C) or (YIL010W and YGR209C)
"""
function parse_grr(s::String)
    gene_list_list = Array{Array{String, 1}, 1}()
    or_genes = split(s, " or ")
    for or_gene in or_genes
        and_genes = split(replace(replace(or_gene, "("=>""), ")"=>""), " and ")
        push!(gene_list_list, and_genes)
    end
    return gene_list_list
end

"""
    save_model(model::CobraTools.Model, file_location::String)

Save model(::CobraTools.Model) at file_location. Infers format from file_location extension.
Supported formats include SBML (.xml), Matlab COBRA models (.mat) and JSON COBRA models (.json).

Note, only the fields contained in model::CobraTools.Model are saved. Make sure that information isn't
lost between reading a model and writing a model (e.g. check gene reaction rules, notes and annotations).
"""
function save_model(model::CobraTools.Model, file_location::String)
    if endswith(file_location, ".json")
        save_json_model(model, file_location)
    elseif endswith(file_location, ".xml")
        @warn "Not implemented!"
    elseif endswith(file_location, ".mat")
        save_matlab_model(model, file_location)
    else
        @error "Model format not supported. The format is inferred from the file extension. Supported formats: *.mat, *.xml, *.json."
    end
end

"""
    save_json_model(model::CobraTools.Model, file_location::String)
"""
function save_json_model(model::CobraTools.Model, file_location::String)
    modeldict = Dict{String, Any}()
    modeldict["id"] = model.id
    modeldict["metabolites"] = model.mets
    modeldict["genes"] = model.genes    
    rxns = []
    for r in model.rxns
        rdict = Dict()
        rdict["id"] = r.id
        rdict["name"] = r.name
        rdict["metabolites"] = Dict{String, Float64}(k.id=>v for (k, v) in r.metabolites)
        rdict["lower_bound"] = r.lb
        rdict["upper_bound"] = r.ub
        rdict["gene_reaction_rule"] = r.grr
        rdict["subsystem"] = r.subsystem
        rdict["notes"] = r.notes
        rdict["annotation"] = r.annotation
        rdict["objective_coefficient"] = r.objective_coefficient
        push!(rxns, rdict)
    end

    modeldict["reactions"] = rxns
    open(file_location, "w") do io
        JSON.print(io, modeldict)
    end
end

"""
    save_matlab_model(model::CobraTools.Model, file_location::String)
"""
function save_matlab_model(model::CobraTools.Model, file_location::String)
    rxnrevs = zeros(Int64, length(model.rxns))
    for i in eachindex(model.rxns)
        if model.rxns[i].lb < 0.0 && model.rxns[i].ub > 0
            rxnrevs[i] = 1 # reversible
        end
    end

    rgm = spzeros(length(model.rxns), length(model.genes)) # stored as a sparse matrix
    for (i, rxn) in enumerate(model.rxns)
        for (j, gene) in enumerate(model.genes)
            if contains(rxn.grr, gene.id)
                rgm[i, j] = 1.0
            end
        end
    end
    
    S, b, ubs, lbs = get_core_model(model)

    mdict = Dict("c" => [r.objective_coefficient for r in model.rxns],
    "mets" => [m.id for m in model.mets],
    "subSystems" => [r.subsystem for r in model.rxns],
    "b" => Array(b),
    "metFormulas" => [m.formula for m in model.mets],
    "rxnGeneMat" => rgm,
    "ub" => Array(ubs),
    "rxnNames" => [r.name for r in model.rxns],
    "description" => model.id,
    "genes" => [g.id for g in model.genes],
    "rev" => rxnrevs,
    "grRules" => [r.grr for r in model.rxns],
    "S" => Array(S),
    "metNames" => [m.name for m in model.mets],
    "lb" => Array(lbs),
    "metCharge" => [m.charge for m in model.mets],
    "rxns" => [r.id for r in model.rxns])

    matwrite(file_location, Dict(Symbol(model.id) => mdict)) 
end

"""
    save_sbml_model(model::CobraTools.Model, file_location::String)
"""
function save_sbml_model(model::CobraTools.Model, file_location::String)
    # To do...
end
