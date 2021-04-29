"""
    check_duplicate_annotations(gene::Gene, genes::Dict{String, Gene}; inspect_annotations...)

Determine if `gene` has any overlapping annotations in `genes`.
The annotations checked are: `inspect_annotations = ["ncbigene", "ncbigi", "refseq_locus_tag", 
"refseq_name", "refseq_synonym", "uniprot"]`.
Return the `id` of the gene with duplicate annotations in `genes`.
If no annotation overlap is found, return `nothing`.
"""
function check_duplicate_annotations(
    check_gene::Gene,
    gs::OrderedDict{String,Gene};
    inspect_annotations = [
        "ncbigene",
        "ncbigi",
        "refseq_locus_tag",
        "refseq_name",
        "refseq_synonym",
        "uniprot",
    ],
)::Union{Nothing,String}
    for (k, gene) in gs
        for anno in inspect_annotations
            if length(
                intersect(
                    get(gene.annotations, anno, ["c1"]),
                    get(check_gene.annotations, anno, "c2"),
                ),
            ) != 0
                return k
            end
        end
    end
    return nothing # no matches found
end
