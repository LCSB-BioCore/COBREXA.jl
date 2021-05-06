
_annotations(m::Metabolite) = m.annotations
_annotations(r::Reaction) = r.annotations
_annotations(g::Gene) = g.annotations

"""
    annotation_index(
        xs::AbstractDict{String};
        annotations = _annotations,
    )::Dict{String,Dict{String,[String]}}

Extract annotations from a dictionary of items `xs` and build an index that
maps annotation "kinds" (e.g. `"PubChem"`) to the mapping from the annotations
(e.g.  `"COMPOUND_12345"`) to item IDs that carry the annotations.

Function `annotations` is used to access the `Annotations` object in the
dictionary values.

This is extremely useful for finding items by annotation data.
"""
function annotation_index(
    xs::AbstractDict{String};
    annotations = _annotations,
)::Dict{String,Dict{String,Set{String}}}
    res = Dict{String,Dict{String,Set{String}}}()
    for (n, ax) in xs
        a = annotations(ax)
        for (k, anns) in a
            if !haskey(res, k)
                res[k] = Dict{String,Set{String}}()
            end
            for v in anns
                if !haskey(res[k], v)
                    res[k][v] = Set([n])
                else
                    push!(res[k][v], n)
                end
            end
        end
    end
    res
end

"""
    ambiguously_identified_items(
        index::Dict{String,Dict{String,[String]}},
    )::Vector{String}

Find items (genes, metabolites, ...) from the annotation index that are
identified non-uniquely by at least one of their annotations.

This often indicates that the items are duplicate or miscategorized.
"""
function ambiguously_identified_items(
    index::Dict{String,Dict{String,Set{String}}},
)::Set{String}
    res = Set{String}()
    for (_, idents) in index
        for (_, items) in idents
            if length(items) > 1
                push!(res, items...)
            end
        end
    end
    res
end
