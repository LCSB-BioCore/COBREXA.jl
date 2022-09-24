"""
$(TYPEDEF)

`StandardModel` is used to store a constraint based metabolic model with
meta-information.  Meta-information is defined as annotation details, which
include gene-reaction-rules, formulas, etc.

This model type seeks to keep as much meta-information as possible, as opposed
to `CoreModel` and `CoreModelCoupled`, which keep the bare neccessities only.
When merging models and keeping meta-information is important, use this as the
model type.  If meta-information is not important, use the more efficient core
model types.  See [`CoreModel`](@ref) and [`CoreModelCoupled`](@ref) for
comparison.

In this model, reactions, metabolites, and genes are stored in ordered
dictionaries indexed by each struct's `id` field.  For example,
`model.reactions["rxn1_id"]` returns a `Reaction` where the field `id` equals
`"rxn1_id"`.  This makes adding and removing reactions efficient.

Note that the stoichiometric matrix (or any other core data, e.g. flux bounds)
is not stored directly as in `CoreModel`.  When this model type is used in
analysis functions, these core data structures are built from scratch each time
an analysis function is called.  This can cause performance issues if you run
many small analysis functions sequentially.  Consider using the core model
types if performance is critical.

See also: [`Reaction`](@ref), [`Metabolite`](@ref), [`Gene`](@ref)

# Example
```
model = load_model(StandardModel, "my_model.json")
keys(model.reactions)
```

# Fields
$(TYPEDFIELDS)
"""
mutable struct StandardModel <: MetabolicModel
    id::String
    reactions::OrderedDict{String,Reaction}
    metabolites::OrderedDict{String,Metabolite}
    genes::OrderedDict{String,Gene}

    StandardModel(
        id = "";
        reactions = OrderedDict{String,Reaction}(),
        metabolites = OrderedDict{String,Metabolite}(),
        genes = OrderedDict{String,Gene}(),
    ) = new(id, reactions, metabolites, genes)
end
