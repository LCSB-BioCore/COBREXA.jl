# Structs and printing

"""
Metabolite
"""
mutable struct Metabolite
    id :: String
    name :: String
    formula :: String
    charge :: Int64
    compartment :: String
    notes :: Dict{String, Array{String, 1}}
    annotation :: Dict{String, Union{Array{String, 1}, String}}
end

function Metabolite()
    id = ""
    name = ""
    formula = ""
    charge = 0
    compartment = ""
    notes = Dict{String, Array{String, 1}}()
    annotation = Dict{String, Union{Array{String, 1}, String}}()
    Metabolite(id, name, formula, charge, compartment, notes, annotation)
end

function Metabolite(id::String)
    name = ""
    formula = ""
    charge = 0
    compartment = ""
    notes = Dict{String, Array{String, 1}}()
    annotation = Dict{String, Union{Array{String, 1}, String}}()
    Metabolite(id, name, formula, charge, compartment, notes, annotation)
end

function Metabolite(d::Dict{String, Any})
    id = ""
    name = ""
    formula = ""
    charge = 0
    compartment = ""
    notes = Dict{String, Array{String, 1}}()
    annotation = Dict{String, Union{Array{String, 1}, String}}()
    for (k, v) in d
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
            notes = Dict{String, Array{String, 1}}(kk=>vv for (kk, vv) in v)
        elseif k == "annotation"
            annotation = Dict{String, Union{Array{String, 1}, String}}()
            for (kk, vv) in v
                if typeof(vv) == String
                    annotation[kk] = vv
                else
                    annotation[kk] = convert(Array{String, 1}, vv)
                end
            end
        else
            cto.verbose && @warn "Unrecognized reaction field: $k"
        end
    end
    
    Metabolite(id, name, formula, charge, compartment, notes, annotation)
end

"""
Reaction
"""
mutable struct Reaction
    id :: String
    name :: String
    metabolites :: Dict{Metabolite, Float64}
    lb :: Float64
    ub :: Float64
    grr :: String
    subsystem :: String
    notes :: Dict{String, Array{String, 1}}
    annotation :: Dict{String, Union{Array{String, 1}, String}} # SBO is a single string term
    objective_coefficient :: Float64
end

function Reaction()
    id = ""
    name = ""
    metabolites = Dict{Metabolite, Float64}()
    lb = -1000.0 
    ub = 1000.0 
    grr = ""
    subsystem = "" 
    notes = Dict{String, Array{String, 1}}() 
    annotation = Dict{String, Union{Array{String, 1}, String}}()
    objective_coefficient = 0.0
    Reaction(id, name, metabolites, lb, ub, grr, subsystem, notes, annotation, objective_coefficient)
end

function Reaction(metabolites::Dict{Metabolite, Float64}, dir::String)
    id = ""
    name = ""
    if dir == "for"
        lb = 0.0 
        ub = 1000.0     
    elseif dir == "rev"
        lb = -1000.0 
        ub = 0.0     
    else
        lb = -1000.0 
        ub = 1000.0     
    end
    grr = ""
    subsystem = "" 
    notes = Dict{String, Array{String, 1}}() 
    annotation = Dict{String, Union{Array{String, 1}, String}}()
    objective_coefficient = 0.0
    Reaction(id, name, metabolites, lb, ub, grr, subsystem, notes, annotation, objective_coefficient)
end

function Reaction(d :: Dict{String, Any}, mets::Array{Metabolite, 1})
    id = ""
    name = ""
    metabolites = Dict{Metabolite, Float64}()
    lb = -1000.0 
    ub = 1000.0 
    grr = ""
    subsystem = "" 
    notes = Dict{String, Array{String, 1}}() 
    annotation = Dict{String, Union{Array{String, 1}, String}}()
    objective_coefficient = 0.0
    for (k, v) in d
        if k == "id"
            id = v
        elseif k == "name"
            name = v
        elseif k == "metabolites"
            metabolites = Dict{Metabolite, Float64}() 
            for (kk, vv) in v
                ind = findfirst(x->x.id == kk, mets)
                isnothing(ind) ? ((cto.verbose && @warn "Metabolite $kk not found in reaction assignment."); continue) : nothing
                metabolites[mets[ind]] = vv
            end 
        elseif k == "lower_bound"
            lb = v
        elseif k == "upper_bound"
            ub = v
        elseif k == "gene_reaction_rule"
            grr = v
        elseif k == "subsystem"
            subsystem = v
        elseif k == "notes"
            notes = Dict{String, Array{String, 1}}(kk=>vv for (kk, vv) in v)
        elseif k == "annotation"
            annotation = Dict{String, Union{Array{String, 1}, String}}()
            for (kk, vv) in v
                if typeof(vv) == String
                    annotation[kk] = vv
                else
                    annotation[kk] = convert(Array{String, 1}, vv)
                end
            end
        elseif k == "objective_coefficient"
            objective_coefficient = v
        else
            cto.verbose && @warn "Unrecognized reaction field: $k"
        end
    end
    Reaction(id, name, metabolites, lb, ub, grr, subsystem, notes, annotation, objective_coefficient)
end


"""
Gene
"""
mutable struct Gene
    id :: String
    name :: String
    notes :: Dict{String, Array{String, 1}}
    annotation :: Dict{String, Union{Array{String, 1}, String}}    
end

function Gene()
    id = ""
    name = ""
    notes = Dict{String, Array{String, 1}}()
    annotation = Dict{String, Union{Array{String, 1}, String}}()    
    Gene(id, name, notes, annotation)
end

function Gene(d)
    id = ""
    name = ""
    notes = Dict{String, Array{String, 1}}()
    annotation = Dict{String, Union{Array{String, 1}, String}}()
    for (k, v) in d
        if k == "id"
            id = v
        elseif k == "name"
            name = v
        elseif k == "notes"
            notes = Dict{String, Array{String, 1}}(kk=>vv for (kk, vv) in v)
        elseif k == "annotation"
            annotation = Dict{String, Union{Array{String, 1}, String}}()
            for (kk, vv) in v
                if typeof(vv) == String
                    annotation[kk] = vv
                else
                    annotation[kk] = convert(Array{String, 1}, vv)
                end
            end
        else
            cto.verbose && @warn "Unrecognized reaction field: $k"
        end
    end
    
    Gene(id, name, notes, annotation)
end

"""
CoreModel fields define the basic information necessary to run analysis tools
"""
struct CoreModel
    S :: SparseMatrixCSC{Float64,Int64} # stoichiometric matrix
    b :: SparseVector{Float64,Int64} # mass balance rhs
    lbs :: Array{Float64, 1} # reaction lower bounds
    ubs :: Array{Float64, 1} # rxn upper bounds
end

"""
CoreModel()

Empty constructor.
"""
function CoreModel()
    CoreModel(sparse(rand(0,0)), sparse(rand(0)), Float64[], Float64[])
end

"""
Complete model with metadata associated with the model e.g. formulas etc.
Contains grrs which should make gene reaction look ups easier
"""
struct Model
    id :: String # model name
    coremodel :: CoreModel # core model
    rxns :: Array{Reaction, 1} # reaction metadata
    mets :: Array{Metabolite, 1}  # metabolite metadata
    genes :: Array{Gene, 1} # gene metadata
    grrs :: Dict{String, Array{Array{String, 1}, 1}} # reaction -> [[gene & gene...] or [gene & gene...] or [gene & gene...]] 
end

"""
Model()

Empty constructor.
"""
function Model()
    Model("blank", CoreModel(), Array{Reaction, 1}(), Array{Metabolite, 1}(), Array{Gene, 1}(), Dict{String, Array{Array{String, 1}, 1}}())
end

"""
index = getindex(model, rxn)

Get the index of rxn in model. Return -1 if not found.
"""
function Base.getindex(model::Model, rxn::Reaction)
    for i in eachindex(model.rxns)
        if model.rxns[i].id == rxn.id
            return i
        end
    end
    return -1
end

"""
index = getindex(model, rxn)

Get the index of metabolite in model. Return -1 if not found.
"""
function Base.getindex(model::Model, met::Metabolite)
    for i in eachindex(model.mets)
        if model.mets[i].id == met.id
            return i
        end
    end
    return -1
end

"""
index = getindex(model, gene)

Get the index of gene in model. Return -1 if not found.
"""
function Base.getindex(model::Model, gene::Gene)
    for i in eachindex(model.genes)
        if model.genes[i].id == gene.id
            return i
        end
    end
    return -1
end

"""
findfirst(rxns, rxnid)

Return the reaction with rxnid or else `nothing`. Typically used: findfirst(model.rxns, rxnid)
"""
function Base.findfirst(rxns::Array{Reaction, 1}, rxnid::String)
    for i in eachindex(rxns)
        if rxns[i].id == rxnid
            return rxns[i]
        end
    end
    return nothing
end

"""
findfirst(mets, metid)

Return the metabolite with metid or else `nothing`. Typically used: findfirst(model.mets, metid)
"""
function Base.findfirst(mets::Array{Metabolite, 1}, metid::String)
    for i in eachindex(mets)
        if mets[i].id == metid
            return mets[i]
        end
    end
    return nothing
end

"""
findfirst(genes, geneid)

Return the gene with geneid or else `nothing`. Typically used: findfirst(model.genes, geneid)
"""
function Base.findfirst(genes::Array{Gene, 1}, geneid::String)
    for i in eachindex(genes)
        if genes[i].id == geneid
            return genes[i]
        end
    end
    return nothing
end


"""
Pretty printing of model.
"""
function Base.show(io::IO, m::Model)
    println(io, "Constraint based model: ", m.id)
    println(io, "Number of reactions: ", length(m.rxns))
    println(io, "Number of metabolites: ", length(m.mets))
    println(io, "Number of genes: ", length(m.genes))
    # println(io, "Objective function: ", )
end

"""
Pretty printing of reaction.
"""
function Base.show(io::IO, r::Reaction)
    if r.ub > 0.0 && r.lb < 0.0
        arrow = " ⟷  "
    elseif r.ub == 0.0 && r.lb < 0.0
        arrow = " ⟵  "
    elseif r.ub > 0.0 && r.lb == 0.0
        arrow = " ⟶  "
    else
        arrow = " →∣←  " # blocked reaction
    end
    substrates = String[]
    products = String[]
    for (k, v) in r.metabolites
        if v < 0.0
            push!(substrates, string(abs(v))*" "*k.id)
        else
            push!(products, string(abs(v))*" "*k.id)
        end
    end
    isempty(substrates) && (substrates = "∅")
    isempty(products) && (products = "∅")
    
    println(io, "Reaction ID: ", r.id)
    println(io, "Reaction name: ", r.name)
    println(io, "Reaction subsystem: ", r.subsystem)
    println(io, join(substrates, " + ")*arrow*join(products, " + "))
    println(io, "Lower bound: ", r.lb)
    println(io, "Upper bound: ", r.ub)
    println(io, "Genes: ", r.grr)
    println(io, "E.C. number: ", join(get(r.annotation, "ec-code", ["N/A"]), " or "))
end

"""
Pretty printing of reactions.
"""
function Base.show(io::IO, ::MIME"text/plain", rs::Array{Reaction, 1})
    println(io, "Reaction set of length: ", length(rs))
end


"""
Pretty printing of metabolite.
"""
function Base.show(io::IO, m::Metabolite)
    println(io, "Metabolite ID: ", m.id)
    println(io, "Metabolite name: ", m.name)
    println(io, "Formula: ", m.formula)
    println(io, "Charge: ", m.charge)
end

"""
Pretty printing of metabolites.
"""
function Base.show(io::IO, ::MIME"text/plain", ms::Array{Metabolite, 1})
    println(io, "Metabolite set of length: ", length(ms))
end

"""
Pretty printing of gene.
"""
function Base.show(io::IO, g::Gene)
    println(io, "Gene ID: ", g.id)
    println(io, "Gene name: ", g.name)
end

"""
Pretty printing of genes.
"""
function Base.show(io::IO, ::MIME"text/plain", gs::Array{Gene, 1})
    println(io, "Gene set of length: ", length(gs))
end
