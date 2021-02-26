"""
Gene struct (mutable).

# Fields
````
id :: String
name :: String
notes :: Dict{String, Array{String, 1}}
annotation :: Dict{String, Union{Array{String, 1}, String}}
````   
"""
mutable struct Gene <: ModelComponent
    id :: String
    name :: String
    notes :: Dict{String, Array{String, 1}}
    annotation :: Dict{String, Union{Array{String, 1}, String}}
end

"""
    Gene()

Empty gene constructor.
"""
function Gene()
    id = ""
    name = ""
    notes = Dict{String, Array{String, 1}}()
    annotation = Dict{String, Union{Array{String, 1}, String}}()    
    Gene(id, name, notes, annotation)
end

"""
    Gene(id::String)

Assign gene with only an id.
"""
function Gene(id::String)
    name = ""
    notes = Dict{String, Array{String, 1}}()
    annotation = Dict{String, Union{Array{String, 1}, String}}()    
    Gene(id, name, notes, annotation)
end

"""
    getindex(genes::Array{Gene, 1}, gene::Gene)

Get the index of a gene in an array of genes. Return -1 if not found.
Typically used `genes[gene] = index`.
"""
function Base.getindex(genes::Array{Gene, 1}, gene::Gene)
    for i in eachindex(genes)
        if genes[i].id == gene.id
            return i
        end
    end
    return -1
end

"""
    findfirst(genes::Array{Gene, 1}, geneid::String)

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
    _is_duplicate(genes::Array{Gene, 1}, gene::Gene)

Check if gene already exists in genes but has another id. 
First check if the ids are the same. 
If not, check if any of the annotations are the same.
"""
function _is_duplicate(genes::Array{Gene, 1}, cgene::Gene)
    for gene in genes
        if gene.id != cgene.id # not same id
            if any([x in get(cgene.annotation, "ncbigene", ["c1"]) for x in get(gene.annotation, "ncbigene", ["c2"])]) ||
                any([x in get(cgene.annotation, "ncbigi", ["c1"]) for x in get(gene.annotation, "ncbigi", ["c2"])]) ||
                any([x in get(cgene.annotation, "refseq_locus_tag", ["c1"]) for x in get(gene.annotation, "refseq_locus_tag", ["c2"])]) ||
                any([x in get(cgene.annotation, "refseq_name", ["c1"]) for x in get(gene.annotation, "refseq_name", ["c2"])]) ||
                any([x in get(cgene.annotation, "refseq_synonym", ["c1"]) for x in get(gene.annotation, "refseq_synonym", ["c2"])]) ||
                any([x in get(cgene.annotation, "uniprot", ["c1"]) for x in get(gene.annotation, "uniprot", ["c2"])])

                return true, genes[gene]
            end
        else
            return true, genes[gene]
        end
    end
    return false, -1
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
function Base.show(io::IO, ::MIME"text/plain", gs::Array{Gene, 1})
    println(io, "Gene set of length: ", length(gs))
end

"""
Pretty printing of grr::Array{Array{Gene,1},1}.
"""
function Base.show(io::IO, ::MIME"text/plain", grr::Array{Array{Gene,1},1})
    grr_strings = String[]
    for gr in grr
        push!(grr_strings, "("*join([g.id for g in gr], " and ")*")")
    end
    println(io, join(grr_strings, " or "))
end
