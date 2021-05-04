"""
Gene struct.

# Fields
````
id :: String
name :: Union{String, Nothing}
notes :: Dict{String, Vector{String}}
annotation :: Dict{String, Union{Vector{String}, String}}
````
"""
mutable struct Gene
    id::String
    name::Maybe{String}
    notes::Notes
    annotations::Annotations

    Gene(id::String = ""; name = nothing, notes = Notes(), annotations = Annotations()) =
        new(id, name, notes, annotations)
end
