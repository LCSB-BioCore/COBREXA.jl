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

    Gene(id::String = ""; name = nothing, notes = Notes(), annotations = Annotations()) =
        new(id, name, notes, annotations)
end
