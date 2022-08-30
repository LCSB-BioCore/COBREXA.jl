
"""
$(TYPEDEF)

A meta-model that represents a model that is serialized on the disk. The
internal model will be loaded on-demand by using any accessor, or by calling
[`precache!`](@ref) directly.

$(TYPEDFIELDS)
"""
mutable struct Serialized{M} <: ModelWrapper where {M<:MetabolicModel}
    m::Maybe{M}
    filename::String

    Serialized{T}(filename::String) where {T} = new{T}(nothing, filename)
    Serialized(model::T, filename::String) where {T<:MetabolicModel} =
        new{T}(model, filename)
end

"""
$(TYPEDSIGNATURES)

Unwrap the serialized model (precaching it transparently).
"""
function unwrap_model(m::Serialized)
    precache!(m)
    m.m
end

"""
$(TYPEDSIGNATURES)

Load the `Serialized` model from disk in case it's not alreadly loaded.
"""
function precache!(model::Serialized)::Nothing
    if isnothing(model.m)
        model.m = deserialize(model.filename)
    end
    nothing
end
