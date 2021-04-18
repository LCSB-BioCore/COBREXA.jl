"""
Gene struct.
This struc is made to contain information used to annotate genes with meta-information.
Here `id` represents the unique ID of the gene, this is the primary access point used to interface with genes in `StandardModel`.
`name` is a shorthand used for readability.
`notes` and `annotation` contain meta-information. Both of these fields are dictionaries that map string IDs to meta-information that
is always stored as a string array EXCEPT for the `"sbo"` field in `annotation` which is just a string. 

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
        id::String = "",
        name::String = "",
        notes = Dict{String,Vector{String}}(),
        annotation = Dict{String,Union{Vector{String},String}}(),
    ) = new(id, name, notes, annotation)
end
