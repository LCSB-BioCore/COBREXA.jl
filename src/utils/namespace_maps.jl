"""
$(TYPEDSIGNATURES)

Returns the first metabolite associated with a reaction, `rid` in `model`.
Optionally, set `use_annotation` to use the first annotation specified by
`use_annotation` instead of the metabolite id. If no annotation is found,
returns the metabolite id.

This function is most useful in conjunction with [`summarize`](@ref).

```
summarize(...; namepace_mapping = reaction_metabolite_map)
summarize(...; namepace_mapping = (model, rid) -> reaction_metabolite_map(model, rid; use_annotation = "bigg.metabolite"))
```
"""
function reaction_metabolite_map(
    model::AbstractMetabolicModel,
    rid::String;
    use_annotation = nothing,
)
    mid = first(keys(reaction_stoichiometry(model, rid)))
    isnothing(use_annotation) && return mid

    anno = first(get(metabolite_annotations(model, mid), use_annotation, [nothing]))
    isnothing(anno) ? mid : anno
end

"""
$(TYPEDSIGNATURES)

Returns the first annotation associated with a reaction, `rid` in `model` from
the namespace denoted by `use_annotation`. If not annotation is found, return
the reaction id.

This function is most useful in conjunction with [`summarize`](@ref).

```
summarize(...; namepace_mapping = reaction_annotation_map)
summarize(...; namepace_mapping = (model, rid) -> reaction_annotation_map(model, rid; use_annotation = "biocyc"))
```
"""
reaction_annotation_map(
    model::AbstractMetabolicModel,
    rid::String;
    use_annotation = "bigg.reaction",
) = first(get(reaction_annotations(model, rid), use_annotation, [rid]))
