"""
Throw an ArgumentError if model.field has any keys in xs::Union{Reaction, Metabolite, Gene}.
"""
function throw_argerror_if_key_found(model, field, xs)
    _ids = collect(keys(getfield(model, field)))
    ids = [x.id for x in xs if x.id in _ids]
    isempty(ids) || throw(ArgumentError("Duplicated $field IDs already present in model: $ids"))
    nothing
end

"""
Throw an ArgumentError if model.field does not have any keys in xs.
"""
function throw_argerror_if_key_missing(model, field, xs::Vector{String})
    _ids = collect(keys(getfield(model, field)))
    ids = [x for x in xs if !(x in _ids)]
    isempty(ids) || throw(ArgumentError("Missing $field IDs in model: $ids"))
    nothing
end
