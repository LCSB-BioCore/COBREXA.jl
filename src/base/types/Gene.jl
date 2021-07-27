"""
Gene struct.

# Fields
````
id :: String
notes :: Dict{String, Vector{String}}
annotation :: Dict{String, Union{Vector{String}, String}}
````
"""
mutable struct Gene
    id::String
    notes::Notes
    annotations::Annotations

    Gene(id::String = ""; notes = Notes(), annotations = Annotations()) =
        new(id, notes, annotations)
end
