
"""
$(TYPEDEF)

A meta-model that represents a model that is serialized on the disk. The
internal model will be loaded on-demand by using any accessor, or by calling
[`precache!`](@ref) directly.

# Fields
$(TYPEDFIELDS)
"""
mutable struct Serialized{M} <: AbstractModelWrapper where {M<:AbstractMetabolicModel}
    m::Maybe{M}
    filename::String

    Serialized{T}(filename::String) where {T} = new{T}(nothing, filename)
    Serialized(model::T, filename::String) where {T<:AbstractMetabolicModel} =
        new{T}(model, filename)
end

"""
$(TYPEDSIGNATURES)
"""
function Accessors.unwrap_model(m::Serialized)
    precache!(m)
    m.m
end

"""
$(TYPEDSIGNATURES)
"""
function Accessors.precache!(model::Serialized)::Nothing
    if isnothing(model.m)
        model.m = deserialize(model.filename)
    end
    nothing
end
