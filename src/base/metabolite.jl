"""
Metabolite struct.

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
    id::String
    name::String
    formula::String
    charge::Int64
    compartment::String
    notes::Dict{String,Array{String,1}}
    annotation::Dict{String,Union{Array{String,1},String}} # everything is a String[] except sbo, which is a String
end

"""
Pretty printing of metabolite::Metabolite.
"""
function Base.show(io::IO, ::MIME"text/plain", m::Metabolite)
    println(
        io,
        "Metabolite ID: ",
        m.id,
        "\n",
        "Metabolite name: ",
        m.name,
        "\n",
        "Formula: ",
        m.formula,
        "\n",
        "Charge: ",
        m.charge,
    )
end

"""
Pretty printing of metabolites::Array{Metabolite, 1}.
"""
function Base.show(io::IO, ::MIME"text/plain", ms::Array{Metabolite,1})
    println(io, "Metabolite set of length: ", length(ms))
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
    notes = Dict{String,Array{String,1}}()
    annotation = Dict{String,Union{Array{String,1},String}}()
    Metabolite(id, name, formula, charge, compartment, notes, annotation)
end

"""
    Metabolite(id::String)

Assigns only the `id` field to a metabolite struct.
"""
function Metabolite(id::String)
    name = ""
    formula = ""
    charge = 0
    compartment = ""
    notes = Dict{String,Array{String,1}}()
    annotation = Dict{String,Union{Array{String,1},String}}()
    Metabolite(id, name, formula, charge, compartment, notes, annotation)
end

"""
    getindex(mets::Array{Metabolite, 1}, met::Metabolite)

Get the index of a `met` in an array of `mets`, based on `id` equality. 
Return -1 if no matches found.
This function overrides the `[]` notation from Base.

Typically used: `index = mets[met]` works. 
"""
function Base.getindex(mets::Array{Metabolite,1}, met::Metabolite)
    for i in eachindex(mets)
        if mets[i].id == met.id
            return i
        end
    end
    return -1
end

"""
    findfirst(mets::Array{Metabolite, 1}, metid::String)

Return the metabolite in `mets` with `metid`, based in `id` field. 
If nothing matches, return `nothing`. 

Typically used: `met = findfirst(model.mets, metid)`.
"""
function Base.findfirst(mets::Array{Metabolite,1}, metid::String)
    for i in eachindex(mets)
        if mets[i].id == metid
            return mets[i]
        end
    end
    return nothing
end

"""
    check_duplicate_annotations(mets::Array{Metabolite, 1}, met::Metabolite)

Determine if a metabolite `met` has overlapping annotations in for some metabolite in `mets`.
If the annotations overlap, then check if they share a compartment to determine if it a a true duplicate.
The annotations checked are: ["kegg.compound", "bigg.metabolite", "chebi", "inchi_key", "sabiork", "hmdb", "seed.compound", "metanetx.chemical", "reactome.compound", "biocyc"].
Return true and the index of the first hit, otherwise false and -1.

See also: [`check_same_formula`](@ref), [`get_atoms`](@ref)
"""
function check_duplicate_annotations(mets::Array{Metabolite,1}, cmet::Metabolite)
    inspect_annotations = [
        "kegg.compound",
        "bigg.metabolite",
        "chebi",
        "inchi_key",
        "sabiork",
        "hmdb",
        "seed.compound",
        "metanetx.chemical",
        "reactome.compound",
        "biocyc",
    ]
    catoms = get_atoms(cmet)
    for met in mets
        for anno in inspect_annotations
            if any([
                x in get(cmet.annotation, anno, ["c1"]) for
                x in get(met.annotation, anno, ["c2"])
            ])
                if met.compartment == cmet.compartment
                    return true, mets[met]
                end
            end
        end
    end
    return false, -1
end

"""
    check_same_formula(mets::Array{Metabolite, 1}, met::Metabolite)

Return all metabolites in `mets` that have the same formula as `met`.
Formula similarity is determined using a hash function, not with [`get_atoms`](@ref).

See also: [`check_duplicate_annotations`](@ref), [`get_atoms`](@ref)
"""
function check_same_formula(mets::Array{Metabolite,1}, met::Metabolite)
    inds = []
    met_formula_hash = hash(met.formula)
    for (i, m) in enumerate(mets)
        if met_formula_hash == hash(m.formula)
            push!(inds, i)
        end
    end
    return mets[inds]
end

"""
    get_atoms(met::Metabolite)

Return a dictionary mapping the elements in a metabolite `met` to their stoichiometric coefficients.

See also: [`check_duplicate_annotations`](@ref), [`check_same_formula`](@ref)
"""
function get_atoms(met::Metabolite)
    atoms = Dict{String,Int64}()
    length(met.formula) == 0 && return atoms
    for m in eachmatch(r"([A-Z]{1})([a-z]?)(\d*)", met.formula)
        element = match(r"([A-Z]{1})([a-z]?)", m.match)
        number = match(r"\d\d*", m.match)
        atoms[string(element.match)] =
            isnothing(number) ? 1 : parse(Int64, string(number.match))
    end
    return atoms
end
