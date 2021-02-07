# Structs and printing

"""
Reaction
"""
mutable struct Reaction
    id :: String
    name :: String
    metabolites :: Dict{String, Float64}
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
    metabolites = Dict{String, Float64}()
    lb = -1000.0 
    ub = 1000.0 
    grr = ""
    subsystem = "" 
    notes = Dict{String, Array{String, 1}}() 
    annotation = Dict{String, Union{Array{String, 1}, String}}()
    objective_coefficient = 0.0
    Reaction(id, name, metabolites, lb, ub, grr, subsystem, notes, annotation, objective_coefficient)
end

function Reaction(id::String, metabolites::Dict{String, Float64}, dir::String)
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

function Reaction(d :: Dict{String, Any})
    id = ""
    name = ""
    metabolites = Dict{String, Float64}()
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
            metabolites = Dict{String, Float64}(kk=>vv for (kk, vv) in v)
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
            @warn "Unrecognized reaction field: $k"
        end
    end
    Reaction(id, name, metabolites, lb, ub, grr, subsystem, notes, annotation, objective_coefficient)
end

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

function Metabolite(id :: String)
    name = ""
    formula = ""
    charge = 0
    compartment = ""
    notes = Dict{String, Array{String, 1}}()
    annotation = Dict{String, Union{Array{String, 1}, String}}()
    Metabolite(id, name, formula, charge, compartment, notes, annotation)
end

function Metabolite(d :: Dict{String, Any})
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
            @warn "Unrecognized reaction field: $k"
        end
    end
    
    Metabolite(id, name, formula, charge, compartment, notes, annotation)
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
            @warn "Unrecognized reaction field: $k"
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

struct CBModel
    model :: JuMP.Model
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
            push!(substrates, string(abs(v))*" "*k)
        else
            push!(products, string(abs(v))*" "*k)
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
function Base.show(io::IO, rs::Array{Reaction, 1})
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
function Base.show(io::IO, ms::Array{Metabolite, 1})
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
function Base.show(io::IO, gs::Array{Gene, 1})
    println(io, "Gene set of length: ", length(gs))
end
