"""
    mutable struct Isozyme

Struct containing isozyme information. Here, `stoichiometry` is a 
dictionary of gene ids to their stoichiometry in the isozyme complex,
and `kcats` is a tuple of the forward and reverse kcats of the isozyme.

# Fields
````
stoichiometry :: Dict{String, Int}
kcats :: Tuple{Float64, Float64}
````
"""
mutable struct Isozyme
    stoichiometry :: Dict{String, Int}
    kcats :: Tuple{Float64, Float64}
end