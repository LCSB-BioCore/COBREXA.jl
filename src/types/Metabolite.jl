"""
Metabolite struct.

# Fields
````
id :: String
name :: String
formula :: String
charge :: Int
compartment :: String
notes :: Dict{String, Vector{String}}
annotation :: Dict{String, Union{Vector{String}, String}}
````
"""
mutable struct Metabolite
    id::String
    name::String
    formula::String
    charge::Int
    compartment::String
    notes::Dict{String,Vector{String}}
    annotation::Dict{String,Vector{String}} # everything is a String[] except sbo, which is a String

    Metabolite(
        id = "";
        name = "",
        formula = "",
        charge = 0,
        compartment = "",
        notes = Dict{String,Vector{String}}(),
        annotation = Dict{String,Vector{String}}(),
    ) = new(id, name, formula, charge, compartment, notes, annotation)
end
