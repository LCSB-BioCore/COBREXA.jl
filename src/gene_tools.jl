"""
Gene struct (mutable)
"""
mutable struct Gene
    id :: String
    name :: String
    notes :: Dict{String, Array{String, 1}}
    annotation :: Dict{String, Union{Array{String, 1}, String}}    
end

"""
gene = Gene()

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
gene = Gene(gene_dict :: Dict{String, Any},)

Assign a gene based on the fields contained in gene_dict.
"""
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
index = getindex(genes::Array{Gene, 1}, gene::Gene)

Get the index of a gene in an array of genes. Return -1 if not found.
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
Pretty printing of gene::Gene.
"""
function Base.show(io::IO, g::Gene)
    println(io, "Gene ID: ", g.id)
    println(io, "Gene name: ", g.name)
end

"""
Pretty printing of genes::Array{Gene, 1}.
"""
function Base.show(io::IO, ::MIME"text/plain", gs::Array{Gene, 1})
    println(io, "Gene set of length: ", length(gs))
end

