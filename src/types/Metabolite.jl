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
    name::Union{String,Nothing}
    formula::Union{String,Nothing}
    charge::Union{Int,Nothing}
    compartment::Union{String,Nothing}
    notes::Dict{String,Vector{String}}
    annotation::Dict{String,Vector{String}} # everything is a String[]

    Metabolite(
        id = "";
        name = nothing,
        formula = nothing,
        charge = nothing,
        compartment = nothing,
        notes = Dict{String,Vector{String}}(),
        annotation = Dict{String,Vector{String}}(),
    ) = new(id, name, formula, charge, compartment, notes, annotation)
end
