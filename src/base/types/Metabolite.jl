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
    name::Maybe{String}
    formula::Maybe{String}
    charge::Maybe{Int}
    compartment::Maybe{String}
    notes::Notes
    annotations::Annotations # everything is a String[]

    Metabolite(
        id = "";
        name = nothing,
        formula = nothing,
        charge = nothing,
        compartment = nothing,
        notes = Dict{String,Vector{String}}(),
        annotations = Dict{String,Vector{String}}(),
    ) = new(id, name, formula, charge, compartment, notes, annotations)
end
