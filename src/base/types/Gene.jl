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
    notes::Maybe{Notes}
    annotations::Maybe{Annotations} # everything is a String[]

    Gene(
        id::String = "";
        name = nothing,
        notes = Dict{String,Vector{String}}(),
        annotations = Dict{String,Union{Vector{String},String}}(),
    ) = new(id, name, notes, annotations)
end
