
"""
    getindex(genes::Vector{Gene}, gene::Gene)

Get the index of a `gene` in an array of `genes` based on `id` field.
Return -1 if no matches found.

Typically used `index = genes[gene]`.
"""
function Base.getindex(genes::Vector{Gene}, gene::Gene)
    for i in eachindex(genes)
        if genes[i].id == gene.id
            return i
        end
    end
    return -1
end

"""
    findfirst(genes::Vector{Gene}, geneid::String)

Return the gene with `geneid` in `genes` or else `nothing`.
Based on matching `id` fields.

Typically used: `gene = findfirst(model.genes, geneid)`.
"""
function Base.findfirst(genes::Vector{Gene}, geneid::String)
    for i in eachindex(genes)
        if genes[i].id == geneid
            return genes[i]
        end
    end
    return nothing
end

"""
    check_duplicate_annotations(genes::Vector{Gene}, gene::Gene)

Determine if a `gene` is has overlapping annotations in `genes`.
The annotations checked are: ["ncbigene", "ncbigi", "refseq_locus_tag", "refseq_name", "refseq_synonym", "uniprot"].
Return true and the index of the first hit, otherwise false and -1.
"""
function check_duplicate_annotations(genes::Vector{Gene}, cgene::Gene)
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

"""
    parsegrr(string_rule, genes::Array{Gene, 1})

Parse a gene reaction rule string `string_rule` into a nested `gene` array `Array{Array{Gene, 1}, 1}`.

Format: (YIL010W and YLR043C) or (YIL010W and YGR209C) where `or` can also be `OR, |, ||` and where `and` can also be `AND, &, &&`.
"""
function _parse_grr(s::String, genes::Array{Gene,1})
    if s == "" || isnothing(s)
        return Array{Array{Gene,1},1}()
    end
    # first get the gene id list in string format
    gene_string_rules = Array{Array{String,1},1}()
    or_genes = split(s, r"\s?(or|OR|(\|\|)|\|)\s?") # separate or terms
    for or_gene in or_genes
        and_genes = split(replace(or_gene, r"\(|\)" => ""), r"\s?(and|AND|(\&\&)|\&)\s?")
        push!(gene_string_rules, and_genes)
    end
    # now map these gene string ids to genes
    grr = Array{Array{Gene,1},1}()
    for gsr in gene_string_rules
        gene_list = Array{Gene,1}()
        for g in gsr
            gene = findfirst(genes, g)
            isnothing(gene) && (@warn "Gene not found..."; continue)
            push!(gene_list, gene)
        end
        push!(grr, gene_list)
    end
    return grr
end

"""
    unparse_grr(grr::Array{Array{Gene, 1}, 1}

Converts a nested `gene` array, `grr`, back into a grr string.
"""
function _unparse_grr(grr::Array{Array{Gene,1},1})
    grr_strings = String[]
    for gr in grr
        push!(grr_strings, "(" * join([g.id for g in gr], " and ") * ")")
    end
    grr_string = join(grr_strings, " or ")
    return grr_string
end
