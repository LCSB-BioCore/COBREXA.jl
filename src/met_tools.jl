"""
Metabolite struct (mutable)
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

"""
metabolite = Metabolite()

Empty metabolite constructor.
"""
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

"""
Metabolite(id::String)

Assigns only the id field to a metabolite struct.
"""
function Metabolite(id::String)
    name = ""
    formula = ""
    charge = 0
    compartment = ""
    notes = Dict{String, Array{String, 1}}()
    annotation = Dict{String, Union{Array{String, 1}, String}}()
    Metabolite(id, name, formula, charge, compartment, notes, annotation)
end

"""
metabolite = Metabolite(field_dict::Dict{String, Any})

Assign a metabolite using fields contained in d.
"""
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
            @warn "Unrecognized reaction field: $k"
        end
    end
    Metabolite(id, name, formula, charge, compartment, notes, annotation)
end

"""
index = getindex(mets::Array{Metabolite, 1}, met::Metabolite)

Get the index of a metabolite in an array of metabolites. Return -1 if not found.
"""
function Base.getindex(mets::Array{Metabolite, 1}, met::Metabolite)
    for i in eachindex(mets)
        if mets[i].id == met.id
            return i
        end
    end
    return -1
end

"""
findfirst(mets::Array{Metabolite, 1}, metid::String)

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
Pretty printing of metabolite::Metabolite.
"""
function Base.show(io::IO, m::Metabolite)
    println(io, "Metabolite ID: ", m.id)
    println(io, "Metabolite name: ", m.name)
    println(io, "Formula: ", m.formula)
    println(io, "Charge: ", m.charge)
end

"""
Pretty printing of metabolites::Array{Metabolite, 1}.
"""
function Base.show(io::IO, ::MIME"text/plain", ms::Array{Metabolite, 1})
    println(io, "Metabolite set of length: ", length(ms))
end

"""
get_atoms(met::Metabolite)

Return a dictionary mapping the elements in a metabolite to their stoichiometric coefficients.
"""
function get_atoms(met::Metabolite)
    atoms = Dict{String, Int64}()
    N = length(met.formula)
    
    N == 0 && return atoms

    caps = findall(x-> isuppercase(x[1]), split(met.formula, ""))
    if length(caps) == 1
        atom, count = split_formula(met.formula)
        atoms[atom] = count
    else
        pend = [caps[2:end].-1;length(met.formula)]
        for (i, j) in zip(caps, pend) 
            atom, count = split_formula(met.formula[i:j])
            atoms[atom] = count
        end
    end
    return atoms
end

"""
atom, count = split_formula(formula)

Split the Atom from the stoichiometric coefficient. 
E.g. C12 => C, 12 or Ca3 => Ca, 3
"""
function split_formula(formula)
    N = length(formula)
    if N > 1
        if islowercase(formula[2][1])
            if N > 2
                atom = string(formula[1:2])
                count = parse(Int64, formula[3:end])
            else
                atom = string(formula[1:2])
                count = 1
            end
        else
            atom = string(formula[1])
            count = parse(Int64, formula[2:end]) 
        end
    else
        atom = string(formula[1])
        count = 1
    end
    return atom, count
end