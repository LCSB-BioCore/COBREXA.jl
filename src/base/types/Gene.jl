"""
Gene struct.

# Fields
````
id :: String
name :: Maybe{String}
notes :: Dict{String, Vector{String}}
annotation :: Dict{String, Union{Vector{String}, String}}
molar_mass :: Maybe{Float64}
````
"""
mutable struct Gene
    id::String
    name::Maybe{String}
    notes::Notes
    annotations::Annotations
    molar_mass::Maybe{Float64}

    Gene(
        id::String = "";
        name = nothing,
        notes = Notes(),
        annotations = Annotations(),
        molar_mass = nothing,
    ) = new(id, name, notes, annotations, molar_mass)
end
