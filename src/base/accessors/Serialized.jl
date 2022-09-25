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
