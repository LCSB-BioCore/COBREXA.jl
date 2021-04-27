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
    name::Union{String, Nothing}
    notes::Dict{String,Vector{String}}
    annotation::Dict{String,Vector{String}} # everything is a String[]

    Gene(
        id::String = "";
        name = nothing,
        notes = Dict{String,Vector{String}}(),
        annotation = Dict{String,Union{Vector{String},String}}(),
    ) = new(id, name, notes, annotation)
end
