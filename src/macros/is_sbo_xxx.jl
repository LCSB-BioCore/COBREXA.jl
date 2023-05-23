
"""
@_is_sbo_fn(anno_id, identifier)

A helper for creating functions like `is_sbo_reaction`, `is_sbo_gene`, `is_sbo_metabolite`.
"""
macro _is_sbo_fn(anno_id, identifiers)

    fname = Symbol(:is_sbo_, anno_id)
    annofunc = Symbol(anno_id, :_annotations)
    accessorfunc = Symbol(anno_id, :s)
    body = quote
        begin
            id in $accessorfunc(model) || return false
            anno = $annofunc(model, id)
            for key in sbo_annotation_keys
                if haskey(anno, key)
                    any(in.($identifiers, Ref(anno[key]))) && return true
                end
            end
            return false
        end
    end

    docstring = """
        $fname(
            model::AbstractMetabolicModel,
            id::String;
            sbo_annotation_keys = ["sbo", "SBO"],
        )

    Check if a $(anno_id) can be identified using the `sbo_annotation_keys`. Returns
    false if no hits or if no keys are found.
    """
    esc(
        Expr(
            :macrocall,
            Symbol("@doc"),
            __source__,
            docstring,
            :(
                $fname(
                    model::AbstractMetabolicModel,
                    id::String;
                    sbo_annotation_keys = ["sbo", "SBO"],
                ) = $body
            ),
        ),
    )
end
