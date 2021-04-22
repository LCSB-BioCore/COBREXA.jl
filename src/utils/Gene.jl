"""
    check_duplicate_annotations(gene::Gene, genes::Dict{String, Gene}, )

Determine if `gene` has any overlapping annotations in `genes`.
The annotations checked are: ["ncbigene", "ncbigi", "refseq_locus_tag", "refseq_name", "refseq_synonym", "uniprot"].
Return true and the `id` of the gene with duplicate annotations in `genes`.
If no annotation overlap is found, return false and "".
"""
function check_duplicate_annotations(
    check_gene::Gene,
    genes::OrderedDict{String,Gene},
)::Tuple{Bool,String}
    inspect_annotations = [
        "ncbigene",
        "ncbigi",
        "refseq_locus_tag",
        "refseq_name",
        "refseq_synonym",
        "uniprot",
    ]
    for (k, gene) in genes
        for anno in inspect_annotations
            if length(
                intersect(
                    get(gene.annotation, anno, ["c1"]),
                    get(check_gene.annotation, anno, "c2"),
                ),
            ) != 0
                return true, k
            end
        end
    end
    return false, ""
end

"""
    parsegrr(string_rule)

Parse a gene reaction rule string `string_rule` into a nested gene reaction rule array `Vector{Vector{String}}`.

Format: (YIL010W and YLR043C) or (YIL010W and YGR209C) where `or` can also be `OR, |, ||` and where `and` can also be `AND, &, &&`.
So `"(YIL010W and YLR043C) or (YIL010W and YGR209C)"`` becomes `[[YIL010W, YLR043C], [YIL010W, YGR209C]]`` as a nested string array.
"""
function _parse_grr(s::String)
    if s == "" || isnothing(s)
        return Vector{Vector{String}}()
    end
    # first get the gene id list in string format
    gene_reaction_rules = Vector{Vector{String}}()
    or_genes = split(s, r"\s?(or|OR|(\|\|)|\|)\s?") # separate or terms
    for or_gene in or_genes
        and_genes = split(replace(or_gene, r"\(|\)" => ""), r"\s?(and|AND|(\&\&)|\&)\s?")
        push!(gene_reaction_rules, and_genes)
    end
    return gene_reaction_rules
end

"""
    unparse_grr(grr::Vector{Vector{Gene}}

Converts a nested gene reaction array  back into a grr string.
So `[[YIL010W, YLR043C], [YIL010W, YGR209C]]`` becomes `"(YIL010W and YLR043C) or (YIL010W and YGR209C)"``.
"""
function _unparse_grr(grr::Vector{Vector{String}})
    grr_strings = String[]
    for gr in grr
        push!(grr_strings, "(" * join([g for g in gr], " and ") * ")")
    end
    grr_string = join(grr_strings, " or ")
    return grr_string
end
