
"""
$(TYPEDSIGNATURES)

Serialize the `model` to file `filename`, returning a [`Serialized`](@ref)
model that can be loaded back transparently by [`precache!`](@ref). The result
does _not_ contain the actual model data that are deferred to the disk; it may
thus be used to save memory, or send the model efficiently to remote workers
within distributed shared-storage environments.

The benefit of using this over "raw" `Serialization.serialize` is that the
resulting `Serialized` model will reload itself automatically with
[`precache!`](@ref) at first usage, which needs to be done manually when using
the `Serialization` package directly.
"""
function serialize_model(
    model::MM,
    filename::String,
)::Serialized{MM} where {MM<:AbstractMetabolicModel}
    open(f -> serialize(f, model), filename, "w")
    Serialized{MM}(filename)
end

"""
$(TYPEDSIGNATURES)

Specialization of [`serialize_model`](@ref) that prevents nested serialization
of already-serialized models.
"""
function serialize_model(model::Serialized, filename::String)
    precache!(model)
    serialize_model(model.m, filename)
end
