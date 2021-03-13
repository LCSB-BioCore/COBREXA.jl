"""
Gene struct.

# Fields
````
id :: String
name :: String
notes :: Dict{String, Array{String, 1}}
annotation :: Dict{String, Union{Array{String, 1}, String}}
````   
"""
mutable struct Gene <: ModelComponent
    id::String
    name::String
    notes::Dict{String,Array{String,1}}
    annotation::Dict{String,Union{Array{String,1},String}} # everything is a String[] except sbo, which is a String
end

"""
Pretty printing of gene::Gene.
"""
function Base.show(io::IO, ::MIME"text/plain", g::Gene)
    println(io, "Gene ID: ", g.id)
    println(io, "Gene name: ", g.name)
end

"""
Pretty printing of genes::Array{Gene, 1}.
"""
function Base.show(io::IO, ::MIME"text/plain", gs::Array{Gene,1})
    println(io, "Gene set of length: ", length(gs))
end

"""
Pretty printing of grr::Array{Array{Gene,1},1}.
"""
function Base.show(io::IO, ::MIME"text/plain", grr::Array{Array{Gene,1},1})
    grr_strings = String[]
    for gr in grr
        push!(grr_strings, "(" * join([g.id for g in gr], " and ") * ")")
    end
    println(io, join(grr_strings, " or "))
end

"""
    Gene()

Empty gene constructor.

See also: [`Gene(::String)`](@ref).
"""
function Gene()
    id = ""
    name = ""
    notes = Dict{String,Array{String,1}}()
    annotation = Dict{String,Union{Array{String,1},String}}()
    Gene(id, name, notes, annotation)
end

"""
    Gene(id::String)

Assign gene with only an `id`.

See also: [`Gene()`](@ref).
"""
function Gene(id::String)
    name = ""
    notes = Dict{String,Array{String,1}}()
    annotation = Dict{String,Union{Array{String,1},String}}()
    Gene(id, name, notes, annotation)
end

"""
    getindex(genes::Array{Gene, 1}, gene::Gene)

Get the index of a `gene` in an array of `genes` based on `id` field. 
Return -1 if no matches found.

Typically used `index = genes[gene]`.
"""
function Base.getindex(genes::Array{Gene,1}, gene::Gene)
    for i in eachindex(genes)
        if genes[i].id == gene.id
            return i
        end
    end
    return -1
end

"""
    findfirst(genes::Array{Gene, 1}, geneid::String)

Return the gene with `geneid` in `genes` or else `nothing`.
Based on matching `id` fields. 

Typically used: `gene = findfirst(model.genes, geneid)`.
"""
function Base.findfirst(genes::Array{Gene,1}, geneid::String)
    for i in eachindex(genes)
        if genes[i].id == geneid
            return genes[i]
        end
    end
    return nothing
end

"""
    check_duplicate_annotations(genes::Array{Gene, 1}, gene::Gene)

Determine if a `gene` is has overlapping annotations in `genes`.
The annotations checked are: ["ncbigene", "ncbigi", "refseq_locus_tag", "refseq_name", "refseq_synonym", "uniprot"].
Return true and the index of the first hit, otherwise false and -1.
"""
function check_duplicate_annotations(genes::Array{Gene,1}, cgene::Gene)
    inspect_annotations = [
        "ncbigene",
        "ncbigi",
        "refseq_locus_tag",
        "refseq_name",
        "refseq_synonym",
        "uniprot",
    ]
    for gene in genes
        for anno in inspect_annotations
            if any([
                x in get(cgene.annotation, anno, ["c1"]) for
                x in get(gene.annotation, anno, ["c2"])
            ])
                return true, genes[gene]
            end
        end
    end
    return false, -1
end
