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
