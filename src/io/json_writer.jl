"""
    save_json_model(model::MetabolicModel, file_name::String)

Save a [`JSONModel`](@ref) in `model` to a JSON file `file_name`.

In case the `model` is not `JSONModel`, it will be converted automatically.
"""
function save_json_model(model::MetabolicModel, file_name::String)
    m =
        typeof(model) == JSONModel ? model :
        begin
            @_io_log @warn "Automatically converting $(typeof(model)) to JSONModel for saving, information may be lost."
            convert(JSONModel, model)
        end

    JSON.print(open(file_name, "w"), m.m)
end

#TODO this needs to be subsumed by the above function with auto-conversion
function _write_model(model::StandardModel, ::Type{JSONFile}, file_location::String)
    modeldict = Dict{String,Any}()
    modeldict["id"] = model.id

    mets = []
    for m in model.metabolites
        mdict = Dict()
        mdict["id"] = m.id
        mdict["name"] = m.name
        mdict["formula"] = m.formula
        mdict["charge"] = m.charge
        mdict["compartment"] = m.compartment
        mdict["notes"] = m.notes
        mdict["annotation"] = m.annotation
        push!(mets, mdict)
    end
    modeldict["metabolites"] = mets

    genes = []
    for g in model.genes
        gdict = Dict()
        gdict["id"] = g.id
        gdict["name"] = g.name
        gdict["notes"] = g.notes
        gdict["annotation"] = g.annotation
        push!(genes, gdict)
    end
    modeldict["genes"] = genes

    rxns = []
    for r in model.reactions
        rdict = Dict()
        rdict["id"] = r.id
        rdict["name"] = r.name
        rdict["metabolites"] = Dict{String,Float64}(k.id => v for (k, v) in r.metabolites)
        rdict["lower_bound"] = r.lb
        rdict["upper_bound"] = r.ub
        rdict["gene_reaction_rule"] = _unparse_grr(r.grr)
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
