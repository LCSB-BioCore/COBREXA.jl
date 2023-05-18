
"""
@_isa_fn(anno_id, identifier)

A helper for creating functions like `isa_reaction`, `isa_gene`, `isa_metabolite`.
"""
macro _isa_fn(anno_id, identifiers)

    fname = Symbol(:isa_, anno_id)
    annofunc = Symbol(anno_id, :_annotations)
    accessorfunc = Symbol(anno_id, :s)
    body = quote
        begin
            id in $accessorfunc(model) || return false
            anno = $annofunc(model, id)
            for key in annotation_keys
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
            annotation_keys = ["sbo", "SBO"],
        )

    Check if a $(anno_id) can be identified using the `annotation_keys`. Returns
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
                    annotation_keys = ["sbo", "SBO"],
                ) = $body
            ),
        ),
    )
end
