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
    notes::Maybe{Notes}
    annotations::Maybe{Annotations} # everything is a String[]

    Metabolite(
        id = "";
        name = nothing,
        formula = nothing,
        charge = nothing,
        compartment = nothing,
        notes = nothing,
        annotations = nothing,
    ) = new(id, name, formula, charge, compartment, notes, annotations)
end
