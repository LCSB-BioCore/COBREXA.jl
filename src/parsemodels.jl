"""
model = readmodel(file_location)

Reads a model file at file_location and returns a constraint based model.
Supported formats include SBML (.xml), Matlab COBRA models (.mat) and JSON COBRA models (.json).
"""
function readmodel(file_location)
    if endswith(file_location, ".json")
        @info "Reading a JSON formatted model..."
        model = reconstructmodeljson(JSON.parsefile(file_location))
        @info "Done reading JSON model."
    elseif endswith(file_location, ".xml")
        @info "Reading an SBML formatted model..."

        @info "Done reading SBML model."
    elseif endswith(file_location, ".mat")
        @info "Reading a Matlab formatted model..."

        @info "Done reading Matlab model."
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
    rxnids = [rxn["id"] for rxn in modeldict["reactions"]]
    ubs = [rxn["upper_bound"] for rxn in modeldict["reactions"]]
    lbs = [rxn["lower_bound"] for rxn in modeldict["reactions"]]
    metids = [met["id"] for met in modeldict["metabolites"]]
    b = spzeros(length(metids))
    grrs = Dict{String,  Array{Array{String, 1}, 1}}()
    S = spzeros(length(metids), length(rxnids))
    
    for (i, rxn) in enumerate(modeldict["reactions"])
        for (met, coeff) in rxn["metabolites"]
            j = findfirst(x -> x == met, metids)
            isnothing(j) ? (@error "S matrix construction error: $met not defined."; continue) : nothing
            S[j, i] = coeff
        end

        isempty(rxn["gene_reaction_rule"]) ? continue : (grrs[rxn["id"]] = parsegrr(rxn["gene_reaction_rule"]))
    end
    
    return Model(id, S, b, lbs, ubs, rxnids, metids, grrs)
end

"""

"""
function reconstructmodelmatlab(file_location)
    mf = MatFile(file_location)
    model_name = variable_names(mf)[1] # assume model name is the only variable
    modeldict = get_variable(mf, model_name)
    close(mf)

    metids = modeldict["mets"]
    rxnids = modeldict["rxns"]
    grrs = Dict{String,  Array{Array{String, 1}, 1}}()
    for (i, rxn) in enumerate(rxnids)
        if !isempty(modeldict["grRules"][i])
            grrs[rxn] = parsegrr(modeldict["grRules"][i])
        end
    end
    S = sparse(modeldict["S"])
    b = sparse(modeldict["b"])
    lbs = modeldict["lb"]
    ubs = modeldict["ub"]
    Model(model_name, S, b, lbs, ubs, rxnids, metids, grrs)
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
    elseif endswith(file_location, ".xml")
        @info "Saving an SBML formatted model..."
    elseif endswith(file_location, ".mat")
        @info "Saving a Matlab formatted model..."
    else
        @error "Model format not supported. The format is inferred from the file extension. Supported formats: *.mat, *.xml, *.json."
    end
end
