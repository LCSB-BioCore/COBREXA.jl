"""
    check_duplicate_annotations(gene::Gene, genes::Dict{String, Gene}; inspect_annotations=_constants.gene_annotation_checks)

Determine if `gene` has any overlapping annotations in `genes`.
The annotations checked are listed in `COBREXA._constants.gene_annotation_checks`.
Return the `id` of the gene with duplicate annotations in `genes`, otherwise `nothing`.
"""
function check_duplicate_annotations(
    check_gene::Gene,
    gs::OrderedDict{String,Gene};
    inspect_annotations = _constants.gene_annotation_checks,
)::Union{Nothing,String}
    for (k, gene) in gs
        if k != check_gene.id
            for anno in inspect_annotations
                if !isempty(intersect(get(gene.annotations, anno, ["c1"]), get(check_gene.annotations, anno, "c2")))
                    return k
                end
            end
        end
    end
    return nothing # no matches found
end
