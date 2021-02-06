"""
model = readmodel(file_location)

Reads a model file at file_location and returns a constraint based model.
Supported formats include SBML (.xml), Matlab COBRA models (.mat) and JSON COBRA models (.json).
"""
function readmodel(file_location)
    if endswith(file_location, ".json")
        @info "Reading a JSON formatted model..."
        try 
            model = reconstructmodeljson(JSON.parsefile(file_location))
            @info "Done reading JSON model."
        catch err
            @error "JSON model reading error.\n$err"
            model = Model()
        end
    elseif endswith(file_location, ".xml")
        @info "Reading an SBML formatted model..."
        try
            model = reconstructmodelsbml(file_location)
            @info "Done reading SBML model."
        catch err
            @error "SBML model reading error.\n$err"
            model = Model()
        end
    elseif endswith(file_location, ".mat")
        @info "Reading a Matlab formatted model..."
       try
            model = reconstructmodelmatlab(file_location)
            @info "Done reading Matlab model." 
       catch err
            @error "Matlab model reading error.\n$err"
            model = Model()
       end
    else
        @error "Model format not supported. The format is inferred from the file extension. Supported formats: *.mat, *.xml, *.json."
        model = Model()
    end
    return model
end

"""
reconstructmodeljson(modeldict)
"""
function reconstructmodeljson(modeldict)
    id = modeldict["id"]

    rxns = Array{Reaction, 1}()
    for rxn in modeldict["reactions"]
        r = Reaction(rxn)
        push!(rxns, r)
    end

    mets = Array{Metabolite, 1}()
    for met in modeldict["metabolites"]
        m = Metabolite(met)
        push!(mets, m)
    end
    
    genes = Array{Gene, 1}()
    for gene in modeldict["genes"]
        g = Gene(gene)
        push!(genes, g)
    end

    ubs = [rxn.ub for rxn in rxns]
    lbs = [rxn.lb for rxn in rxns]
    
    b = spzeros(length(mets))
    S = spzeros(length(mets), length(rxns))
    grrs = Dict{String,  Array{Array{String, 1}, 1}}()

    metids = [met.id for met in mets] # need indices for S matrix construction
    for (i, rxn) in enumerate(rxns) # column
        for (met, coeff) in rxn.metabolites
            j = findfirst(x -> x == met, metids) # row
            isnothing(j) ? (@error "S matrix construction error: $met not defined."; continue) : nothing
            S[j, i] = coeff
        end

        isempty(rxn.grr) ? continue : (grrs[rxn.id] = parsegrr(rxn.grr))
    end
    
    return Model(id, CoreModel(S, b, lbs, ubs), rxns, mets, genes, grrs)
end

"""
reconstructmodelmatlab(file_location)
"""
function reconstructmodelmatlab(file_location)
    mf = MatFile(file_location)
    model_name = variable_names(mf)[1] # assume model name is the only variable
    modeldict = get_variable(mf, model_name)
    close(mf)

    model_id = haskey(modeldict, "description") ? modeldict["description"] : model_name
    
    S = sparse(modeldict["S"])
    b = sparse(modeldict["b"])
    lbs = modeldict["lb"]
    ubs = modeldict["ub"]

    mets = Array{Metabolite, 1}()
    for i in eachindex(modeldict["mets"])
        id = haskey(modeldict, "mets") ? modeldict["mets"][i] : ""
        name = haskey(modeldict, "metNames") ? modeldict["metNames"][i] : ""
        formula = haskey(modeldict, "metFormulas") ? modeldict["metFormulas"][i] : "" 
        charge = haskey(modeldict, "metCharge") ? modeldict["metCharge"][i] : 0

        # these fields likely don't exist in the matlab model
        compartment = haskey(modeldict, "compartment") ? modeldict["compartment"][i] : ""
        notes = haskey(modeldict, "notes") ? modeldict["notes"][i] : Dict{String, Array{String, 1}}()
        annotation = haskey(modeldict, "annotation") ? modeldict["annotation"][i] : Dict{String, Union{Array{String, 1}, String}}()

        push!(mets, Metabolite(id, name, formula, charge, compartment, notes, annotation))
    end

    rxns = Array{Reaction, 1}()
    for i in eachindex(modeldict["rxns"])
        id = haskey(modeldict, "rxns") ? modeldict["rxns"][i] : ""
        name = haskey(modeldict, "rxnNames") ? modeldict["rxnNames"][i] : ""
        metinds = findall(x -> x .!= 0.0, modeldict["S"][:, i])
        metabolites = Dict{String, Float64}(mets[j].id=>modeldict["S"][j, i] for j in metinds)
        lb = haskey(modeldict, "lb") ? modeldict["lb"][i] : -1000.0 # reversible by default
        ub = haskey(modeldict, "ub") ? modeldict["ub"][i] : 1000.0 # reversible by default
        grr = haskey(modeldict, "grRules") ? modeldict["grRules"][i] : ""
        subsystem = modeldict["subSystems"][i]
        objective_coefficient = haskey(modeldict, "c") ? modeldict["c"][i] : 0.0

        # these fields likely don't exist in the matlab model
        notes = haskey(modeldict, "notes") ? modeldict["notes"][i] : Dict{String, Array{String, 1}}()
        annotation = haskey(modeldict, "annotation") ? modeldict["annotation"][i] : Dict{String, Union{Array{String, 1}, String}}() 
        
        push!(rxns, Reaction(id, name, metabolites, lb, ub, grr, subsystem, notes, annotation, objective_coefficient))
    end

    genes = Array{Gene, 1}()
    for i in eachindex(modeldict["genes"])
        id = haskey(modeldict, "genes") ? modeldict["genes"][i] : ""
        
        # these fields likely don't exist in the matlab model
        name = haskey(modeldict, "geneNames") ? modeldict["geneNames"][i] : ""
        notes = haskey(modeldict, "geneNotes") ? modeldict["geneNotes"][i] : Dict{String, Array{String, 1}}()
        annotation = haskey(modeldict, "geneAnnotations") ? modeldict["geneAnnotations"][i] : Dict{String, Union{Array{String, 1}, String}}()
        
        push!(genes, Gene(id, name, notes, annotation))
    end

    grrs = Dict{String,  Array{Array{String, 1}, 1}}()
    for (i, rxn) in enumerate(rxns)
        if !isempty(modeldict["grRules"][i])
            grrs[rxn.id] = parsegrr(modeldict["grRules"][i])
        end
    end
    
    return Model(model_id, CoreModel(S, b, lbs, ubs), rxns, mets, genes, grrs)
end

function reconstructmodelsbml(file_location)
    @pyimport libsbml
    reader = libsbml.SBMLReader()
    sbmlmodel = reader[:readSBML](file_location)
    
    mets = Array{Metabolite, 1}()
    # for i in eachindex(modeldict["mets"])
    #     id = haskey(modeldict, "mets") ? modeldict["mets"][i] : ""
    #     name = haskey(modeldict, "metNames") ? modeldict["metNames"][i] : ""
    #     formula = haskey(modeldict, "metFormulas") ? modeldict["metFormulas"][i] : "" 
    #     charge = haskey(modeldict, "metCharge") ? modeldict["metCharge"][i] : 0

    #     # these fields likely don't exist in the matlab model
    #     compartment = haskey(modeldict, "compartment") ? modeldict["compartment"][i] : ""
    #     notes = haskey(modeldict, "notes") ? modeldict["notes"][i] : Dict{String, Array{String, 1}}()
    #     annotation = haskey(modeldict, "annotation") ? modeldict["annotation"][i] : Dict{String, Union{Array{String, 1}, String}}()

    #     push!(mets, Metabolite(id, name, formula, charge, compartment, notes, annotation))
    # end

    rxns = Array{Reaction, 1}()
    # for i in eachindex(modeldict["rxns"])
    #     id = haskey(modeldict, "rxns") ? modeldict["rxns"][i] : ""
    #     name = haskey(modeldict, "rxnNames") ? modeldict["rxnNames"][i] : ""
    #     metinds = findall(x -> x .!= 0.0, modeldict["S"][:, i])
    #     metabolites = Dict{String, Float64}(mets[j].id=>modeldict["S"][j, i] for j in metinds)
    #     lb = haskey(modeldict, "lb") ? modeldict["lb"][i] : -1000.0 # reversible by default
    #     ub = haskey(modeldict, "ub") ? modeldict["ub"][i] : 1000.0 # reversible by default
    #     grr = haskey(modeldict, "grRules") ? modeldict["grRules"][i] : ""
    #     subsystem = modeldict["subSystems"][i]
    #     objective_coefficient = haskey(modeldict, "c") ? modeldict["c"][i] : 0.0

    #     # these fields likely don't exist in the matlab model
    #     notes = haskey(modeldict, "notes") ? modeldict["notes"][i] : Dict{String, Array{String, 1}}()
    #     annotation = haskey(modeldict, "annotation") ? modeldict["annotation"][i] : Dict{String, Union{Array{String, 1}, String}}() 
        
    #     push!(rxns, Reaction(id, name, metabolites, lb, ub, grr, subsystem, notes, annotation, objective_coefficient))
    # end

    genes = Array{Gene, 1}()
    # for i in eachindex(modeldict["genes"])
    #     id = haskey(modeldict, "genes") ? modeldict["genes"][i] : ""
        
    #     # these fields likely don't exist in the matlab model
    #     name = haskey(modeldict, "geneNames") ? modeldict["geneNames"][i] : ""
    #     notes = haskey(modeldict, "geneNotes") ? modeldict["geneNotes"][i] : Dict{String, Array{String, 1}}()
    #     annotation = haskey(modeldict, "geneAnnotations") ? modeldict["geneAnnotations"][i] : Dict{String, Union{Array{String, 1}, String}}()
        
    #     push!(genes, Gene(id, name, notes, annotation))
    # end

    # grrs = Dict{String,  Array{Array{String, 1}, 1}}()
    # for (i, rxn) in enumerate(rxns)
    #     if !isempty(modeldict["grRules"][i])
    #         grrs[rxn.id] = parsegrr(modeldict["grRules"][i])
    #     end
    # end

    return Model()
    # return Model(model_id, CoreModel(S, b, lbs, ubs), rxns, mets, genes, grrs)
end

"""
parsegrr(string_rule)

Format: (YIL010W and YLR043C) or (YIL010W and YGR209C)
"""
function parsegrr(s :: String)
    gene_list_list = Array{Array{String, 1}, 1}()
    or_genes = split(s, " or ")
    for or_gene in or_genes
        and_genes = split(replace(replace(or_gene, "("=>""), ")"=>""), " and ")
        push!(gene_list_list, and_genes)
    end
    return gene_list_list
end

"""
savemodel(model, file_location)

Save model at location file_location. Infers format from file_location extension.
Supported formats include SBML (.xml), Matlab COBRA models (.mat) and JSON COBRA models (.json).
"""
function savemodel(model :: Model, file_location :: String)
    if endswith(file_location, ".json")
        @info "Saving a JSON formatted model..."
        @info "Done saving JSON model."
    elseif endswith(file_location, ".xml")
        @info "Saving an SBML formatted model..."
        @info "Done saving SBML model."
    elseif endswith(file_location, ".mat")
        @info "Saving a Matlab formatted model..."
        @info "Done saving Matlab model."
    else
        @error "Model format not supported. The format is inferred from the file extension. Supported formats: *.mat, *.xml, *.json."
    end
end


function savejsonmodel()

end

function savematlabmodel()

end

function savesbmlmodel()

end

