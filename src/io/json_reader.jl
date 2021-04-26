function _read_model(filename::String, ::Type{JSONFile}, ::Type{JSONModel})
    return JSONModel(JSON.parsefile(filename))
end

function _read_model(file_location::String, ::Type{JSONFile}, ::Type{StandardModel})
    #TODO this should be subsumed by the above loader + abstract model + StdModel converter
    modeldict = JSON.parsefile(file_location)

    modelid = modeldict["id"]

    mets = OrderedDict{String,Metabolite}()
    for i in modeldict["metabolites"]
        met = Metabolite()

        for (k, v) in i
            if k == "id"
                met.id = v
            elseif k == "name"
                met.name = v
            elseif k == "formula"
                met.formula = v
            elseif k == "charge"
                met.charge = v
            elseif k == "compartment"
                met.compartment = v
            elseif k == "notes"
                met.notes = Dict{String,Vector{String}}(kk => vv for (kk, vv) in v)
            elseif k == "annotation"
                for (kk, vv) in v
                    if typeof(vv) == String
                        met.annotation[kk] = [vv]
                    else
                        met.annotation[kk] = convert(Vector{String}, vv)
                    end
                end
            else
                @warn "Unrecognized metabolite field: $k"
            end
        end
        mets[met.id] = met
    end

    genes = OrderedDict{String,Gene}()
    for i in modeldict["genes"]
        gene = Gene()
        for (k, v) in i
            if k == "id"
                gene.id = v
            elseif k == "name"
                gene.name = v
            elseif k == "notes"
                for (kk, vv) in v
                    gene.notes[kk] = vv
                end
            elseif k == "annotation"
                for (kk, vv) in v
                    if typeof(vv) == String
                        gene.annotation[kk] = [vv]
                    else
                        gene.annotation[kk] = convert(Vector{String}, vv)
                    end
                end
            else
                @warn "Unrecognized gene field: $k"
            end
        end
        genes[gene.id] = gene
    end

    rxns = OrderedDict{String,Reaction}()
    for i in modeldict["reactions"]
        rxn = Reaction()
        for (k, v) in i
            if k == "id"
                rxn.id = v
            elseif k == "name"
                rxn.name = v
            elseif k == "metabolites"
                for (kk, vv) in v
                    rxn.metabolites[mets[kk].id] = vv
                end
            elseif k == "lower_bound"
                rxn.lb = v
            elseif k == "upper_bound"
                rxn.ub = v
            elseif k == "gene_reaction_rule"
                rxn.grr = _parse_grr(v)
            elseif k == "subsystem"
                rxn.subsystem = v
            elseif k == "notes"
                for (kk, vv) in v
                    rxn.notes[kk] = vv
                end
            elseif k == "annotation"
                for (kk, vv) in v
                    if typeof(vv) == String
                        rxn.annotation[kk] = [vv]
                    else
                        rxn.annotation[kk] = convert(Vector{String}, vv)
                    end
                end
            elseif k == "objective_coefficient"
                rxn.objective_coefficient = v
            else
                @warn "Unrecognized reaction field: $k"
            end
        end
        rxns[rxn.id] = rxn
    end

    return StandardModel(modelid; reactions = rxns, metabolites = mets, genes = genes)
end
