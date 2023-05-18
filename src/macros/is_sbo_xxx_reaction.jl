
"""
@_is_sbo_reaction_fn(anno_id, identifier)

A helper for creating functions like `is_sbo_exchange_reaction`.
"""
macro _is_sbo_reaction_fn(anno_id, identifiers)

    fname = Symbol(:is_sbo_, anno_id, :_reaction)
    grammar = any(startswith.(anno_id, ["a", "e", "i", "o", "u"])) ? "an" : "a"

    body = quote
        begin
            reaction_id in reactions(model) || return false
            anno = reaction_annotations(model, reaction_id)
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
            reaction_id::String;
            sbo_annotation_keys = ["sbo", "SBO"],
        )

    Check if a reaction is annotated as $(grammar) SBO annotated $(anno_id)
    reaction. In the reaction annotations, use the keys in `sbo_annotation_keys`
    to look for entries. Returns false if no hits or if no keys are found.
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
                    sbo_annotation_keys = ["sbo", "SBO"],
                ) = $body
            ),
        ),
    )
end
