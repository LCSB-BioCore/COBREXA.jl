"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
mutable struct Gene
    id::String
    name::Maybe{String}
    notes::Notes
    annotations::Annotations
end

"""
$(TYPEDSIGNATURES)

A convenient constructor for a `Gene`.
"""
Gene(id = ""; name = nothing, notes = Notes(), annotations = Annotations()) =
    Gene(String(id), name, notes, annotations)
