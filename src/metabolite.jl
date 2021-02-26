"""
Metabolite struct (mutable).

# Fields
````
id :: String
name :: String
formula :: String
charge :: Int64
compartment :: String
notes :: Dict{String, Array{String, 1}}
annotation :: Dict{String, Union{Array{String, 1}, String}}
````
"""
mutable struct Metabolite <: ModelComponent
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
    getindex(mets::Array{Metabolite, 1}, met::Metabolite)

Get the index of a metabolite in an array of metabolites. Return -1 if not found.
This function overrides the [] notation from base, hence `mets[met] = index` works. 
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
    _is_duplicate(mets::Array{Metabolite, 1}, met::Metabolite)

Check if met already exists in mets but has another id. 
First check if the id and formulas are the same. 
If not, check if the charges are the same.
If not, check if any of the annotations are the same.
"""
function _is_duplicate(mets::Array{Metabolite, 1}, cmet::Metabolite)
    catoms = get_atoms(cmet)
    for met in mets
        if met.id != cmet.id
            matoms = get_atoms(met)
            if all([matoms[k]==v for (k, v) in catoms]) # all the atoms are the same
                if cmet.charge == met.charge # if charges the same
                    if any([x in get(cmet.annotation, "kegg.compound", ["c1"]) for x in get(met.annotation, "kegg.compound", ["c2"])]) ||
                        any([x in get(cmet.annotation, "bigg.metabolite", ["c1"]) for x in get(met.annotation, "bigg.metabolite", ["c2"])]) ||
                        any([x in get(cmet.annotation, "chebi", ["c1"]) for x in get(met.annotation, "chebi", ["c2"])]) ||
                        any([x in get(cmet.annotation, "inchi_key", ["c1"]) for x in get(met.annotation, "inchi_key", ["c2"])]) ||
                        any([x in get(cmet.annotation, "sabiork", ["c1"]) for x in get(met.annotation, "sabiork", ["c2"])]) ||
                        any([x in get(cmet.annotation, "hmdb", ["c1"]) for x in get(met.annotation, "hmdb", ["c2"])]) ||
                        any([x in get(cmet.annotation, "seed.compound", ["c1"]) for x in get(met.annotation, "seed.compound", ["c2"])]) ||
                        any([x in get(cmet.annotation, "metanetx.chemical", ["c1"]) for x in get(met.annotation, "metanetx.chemical", ["c2"])]) 
                        
                        return true, mets[met]
                    end
                end
            end
        else
            return true, mets[met]
        end
    end
    return false, -1
end

"""
Pretty printing of metabolite::Metabolite.
"""
function Base.show(io::IO, ::MIME"text/plain", m::Metabolite)
    println(io, "Metabolite ID: ", m.id, "\n",
                "Metabolite name: ", m.name, "\n",
                "Formula: ", m.formula, "\n",
                "Charge: ", m.charge)
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