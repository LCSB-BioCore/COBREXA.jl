"""
    read_json_model(modeldict::String)
"""
function read_json_model(file_location::String)
    modeldict = JSON.parsefile(file_location)

    modelid = modeldict["id"]

    mets = Metabolite[]
    for met in modeldict["metabolites"]
        id = ""
        name = ""
        formula = ""
        charge = 0
        compartment = ""
        notes = Dict{String,Array{String,1}}()
        annotation = Dict{String,Union{Array{String,1},String}}()
        for (k, v) in met
            if k == "id"
                id = v
            elseif k == "name"
                name = v
            elseif k == "formula"
                formula = v
            elseif k == "charge"
                charge = v
            elseif k == "compartment"
                compartment = v
            elseif k == "notes"
                notes = Dict{String,Array{String,1}}(kk => vv for (kk, vv) in v)
            elseif k == "annotation"
                for (kk, vv) in v
                    if typeof(vv) == String
                        annotation[kk] = vv
                    else
                        annotation[kk] = convert(Array{String,1}, vv)
                    end
                end
            else
                @warn "Unrecognized reaction field: $k"
            end
        end
        push!(mets, Metabolite(id, name, formula, charge, compartment, notes, annotation))
    end

    genes = Gene[]
    for gene in modeldict["genes"]
        id = ""
        name = ""
        notes = Dict{String,Array{String,1}}()
        annotation = Dict{String,Union{Array{String,1},String}}()
        for (k, v) in gene
            if k == "id"
                id = v
            elseif k == "name"
                name = v
            elseif k == "notes"
                notes = Dict{String,Array{String,1}}(kk => vv for (kk, vv) in v)
            elseif k == "annotation"
                for (kk, vv) in v
                    if typeof(vv) == String
                        annotation[kk] = vv
                    else
                        annotation[kk] = convert(Array{String,1}, vv)
                    end
                end
            else
                @warn "Unrecognized reaction field: $k"
            end
        end
        push!(genes, Gene(id, name, notes, annotation))
    end

    rxns = Reaction[]
    for rxn in modeldict["reactions"]
        id = ""
        name = ""
        metabolites = Dict{Metabolite,Float64}()
        lb = -1000.0
        ub = 1000.0
        grr = Array{Array{Gene,1},1}()
        subsystem = ""
        notes = Dict{String,Array{String,1}}()
        annotation = Dict{String,Union{Array{String,1},String}}()
        objective_coefficient = 0.0
        for (k, v) in rxn
            if k == "id"
                id = v
            elseif k == "name"
                name = v
            elseif k == "metabolites"
                metabolites = Dict{Metabolite,Float64}()
                for (kk, vv) in v
                    ind = findfirst(x -> x.id == kk, mets)
                    isnothing(ind) ?
                    (@warn "Metabolite $kk not found in reaction assignment."; continue) :
                    nothing
                    metabolites[mets[ind]] = vv
                end
            elseif k == "lower_bound"
                lb = v
            elseif k == "upper_bound"
                ub = v
            elseif k == "gene_reaction_rule"
                grr = parse_grr(v, genes)
            elseif k == "subsystem"
                subsystem = v
            elseif k == "notes"
                notes = Dict{String,Array{String,1}}(kk => vv for (kk, vv) in v)
            elseif k == "annotation"
                for (kk, vv) in v
                    if typeof(vv) == String
                        annotation[kk] = vv
                    else
                        annotation[kk] = convert(Array{String,1}, vv)
                    end
                end
            elseif k == "objective_coefficient"
                objective_coefficient = v
            else
                @warn "Unrecognized reaction field: $k"
            end
        end
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

    return StandardModel(modelid, rxns, mets, genes)
end
