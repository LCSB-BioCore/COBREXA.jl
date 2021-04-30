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
    annotations::Annotations # everything is a String[]
    reactions::Set{String}

    Gene(
        id::String = "";
        name = nothing,
        notes = Notes(),
        annotations = Annotations(),
        reactions = Set{String}(),
    ) = new(id, name, notes, annotations, reactions)
end
