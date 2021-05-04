"""
Gene struct.

# Fields
````
id :: String
name :: Union{String, Nothing}
notes :: Dict{String, Vector{String}}
annotation :: Dict{String, Union{Vector{String}, String}}
associated_reactions :: Set{String}
````
"""
mutable struct Gene
    id::String
    name::Maybe{String}
    notes::Notes
    annotations::Annotations
    associated_reactions::Set{String}

    Gene(
        id::String = "";
        name = nothing,
        notes = Notes(),
        annotations = Annotations(),
        associated_reactions = Set{String}(),
    ) = new(id, name, notes, annotations, associated_reactions)
end
