
"""
    serialize_model(model::MM, filename::String)::Serialized{MM} where {MM<:MetabolicModel}

Serialize the `model` to file `filename`, returning a [`Serialized`](@ref)
model that is able to load itself back automatically upon precaching by
[`precache!`](@ref).
"""
function serialize_model(
    model::MM,
    filename::String,
)::Serialized{MM} where {MM<:MetabolicModel}
    open(f -> serialize(f, model), filename, "w")
    Serialized{MM}(nothing, filename)
end

"""
    serialize_model(model::Serialized, filename::String)::Serialized

Specialization of [`serialize_model`](@ref) that prevents nested serialization
of already-serialized models.
"""
function serialize_model(model::Serialized, filename::String)
    precache!(model)
    serialize_model(model.m, filename)
end
