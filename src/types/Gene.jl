"""
Gene struct.

# Fields
````
id :: String
name :: String
notes :: Dict{String, Vector{String}}
annotation :: Dict{String, Union{Vector{String}, String}}
````
"""
mutable struct Gene
    id::String
    name::String
    notes::Dict{String,Vector{String}}
    annotation::Dict{String,Union{Vector{String},String}} # everything is a String[] except sbo, which is a String

    Gene(
        id::String = "";
        name::String = "",
        notes = Dict{String,Vector{String}}(),
        annotation = Dict{String,Union{Vector{String},String}}(),
    ) = new(id, name, notes, annotation)
end
