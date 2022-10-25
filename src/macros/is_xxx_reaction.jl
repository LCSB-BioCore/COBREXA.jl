
"""
@_is_reaction_fn(anno_id, identifier)

A helper for creating functions like `is_exchange_reaction`.
"""
macro _is_reaction_fn(anno_id, identifiers)

    fname = Symbol(:is_, anno_id, :_reaction)
    grammar = any(startswith.(anno_id, ["a", "e", "i", "o", "u"])) ? "an" : "a"

    body = quote
        begin
            anno = reaction_annotations(model, reaction_id)
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
            reaction_id::String;
            annotation_keys = ["sbo", "SBO"],
        )

    Check if a reaction is annotated as $(grammar) $(anno_id) reaction. Uses
    `$identifiers` internally, which includes SBO identifiers. In
    the reaction annotations, use the keys in `annotation_keys` to look for entries.
    Returns false if no hits or if no keys are found.
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
                    reaction_id::String;
                    annotation_keys = ["sbo", "SBO"],
                ) = $body
            ),
        ),
    )
end
