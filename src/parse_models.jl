"""
model = readmodel(file_location)

Reads a model file at file_location and returns a constraint based model.
Supported formats include SBML (.xml), Matlab COBRA models (.mat) and JSON COBRA models (.json).

Note, when importing JSON models only reactions, metabolites, genes and id are used. 

Note, when importing Matlab models only rxns, metCharge, lb, metNames, S, grRules,
genes, description, rxnNames, ub, metFormulas, b, subSystems, mets and c are used.

Note, SBML is not implemented yet.
"""
function readmodel(file_location)
    if endswith(file_location, ".json")
        try 
            model = reconstructmodeljson(JSON.parsefile(file_location))
        catch err
            @error "JSON model reading error.\n$err"
            model = Model()
        end
    elseif endswith(file_location, ".xml")
        try
            model = reconstructmodelsbml(file_location)
            cto.verbose && @warn "Not implemented!"
        catch err
            @error "SBML model reading error.\n$err"
            model = Model()
        end
    elseif endswith(file_location, ".mat")
       try
            model = reconstructmodelmatlab(file_location)
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
        isempty(rxn.grr) ? continue : (grrs[rxn.id] = parsegrr(rxn.grr))
    end
    
    return Model(id, rxns, mets, genes, grrs)
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
    
    mets = Metabolite[]
    for i in eachindex(modeldict["mets"])
        id = haskey(modeldict, "mets") ? modeldict["mets"][i] : ""
        name = haskey(modeldict, "metNames") ? modeldict["metNames"][i] : ""
        formula = haskey(modeldict, "metFormulas") ? modeldict["metFormulas"][i] : "" 
        charge = haskey(modeldict, "metCharge") ? modeldict["metCharge"][i] : 0

        # these fields likely don't exist in the matlab model
        compartment = haskey(modeldict, "compartment") ? modeldict["compartment"][i] : ""
        notes = haskey(modeldict, "notes") ? modeldict["notes"][i] : Dict{String, Array{String, 1}}()
        annotation = haskey(modeldict, "annotation") ? modeldict["annotation"][i] : Dict{String, Union{Array{String, 1}, String}}()

        push!(mets, Metabolite(id, name, formula, charge, compartment, notes, annotation, 1e-3)) # concentration 1 mM
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
        subsystem = modeldict["subSystems"][i]
        objective_coefficient = haskey(modeldict, "c") ? modeldict["c"][i] : 0.0

        # these fields likely don't exist in the matlab model
        notes = haskey(modeldict, "notes") ? modeldict["notes"][i] : Dict{String, Array{String, 1}}()
        annotation = haskey(modeldict, "annotation") ? modeldict["annotation"][i] : Dict{String, Union{Array{String, 1}, String}}() 
        
        push!(rxns, Reaction(id, name, metabolites, lb, ub, grr, subsystem, notes, annotation, objective_coefficient))
    end

    genes = Gene[]
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
    
    return Model(model_id, rxns, mets, genes, grrs)
end

function reconstructmodelsbml(file_location)
    # @pyimport libsbml
    # reader = libsbml.SBMLReader()
    # sbmldoc = reader[:readSBML](file_location)
    # sbmlmodel = sbmldoc[:getModel]() # Get the model
    # model_id = sbmlmodel[:getId]()
    # met = sbmlmodel[:getListOfSpecies]()[1]
    # id = met[:getId]()
    # name = met[:getName]()
    # charge = met[:getCharge]()
    # formula = ???
    return Model()
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

Note, when exporting JSON models only reactions, metabolites, genes and id are written. 

Note, when exporting Matlab models only rxns, metCharge, lb, metNames, S, grRules,
genes, description, rxnNames, ub, metFormulas, b, subSystems, mets, rev, rxnGeneMat
and c are written to.

Note, SBML is not implemented yet.
"""
function savemodel(model :: Model, file_location :: String)
    if endswith(file_location, ".json")
        savejsonmodel(model, file_location)
    elseif endswith(file_location, ".xml")
        cto.verbose && @warn "Not implemented!"
    elseif endswith(file_location, ".mat")
        savematlabmodel(model, file_location)
    else
        @error "Model format not supported. The format is inferred from the file extension. Supported formats: *.mat, *.xml, *.json."
    end
end


function savejsonmodel(model :: Model, file_location :: String)
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

function savematlabmodel(model :: Model, file_location :: String)
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
    
    coremodel = CoreModel(model::Model)

    mdict = Dict("c" => [r.objective_coefficient for r in model.rxns],
    "mets" => [m.id for m in model.mets],
    "subSystems" => [r.subsystem for r in model.rxns],
    "b" => Array(coremodel.b),
    "metFormulas" => [m.formula for m in model.mets],
    "rxnGeneMat" => rgm,
    "ub" => Array(coremodel.ubs),
    "rxnNames" => [r.name for r in model.rxns],
    "description" => model.id,
    "genes" => [g.id for g in model.genes],
    "rev" => rxnrevs,
    "grRules" => [r.grr for r in model.rxns],
    "S" => Array(coremodel.S),
    "metNames" => [m.name for m in model.mets],
    "lb" => Array(coremodel.lbs),
    "metCharge" => [m.charge for m in model.mets],
    "rxns" => [r.id for r in model.rxns])

    write_matfile(file_location; Dict(Symbol(model.id) => mdict)...) 
end

function savesbmlmodel()
    # To do...
end

