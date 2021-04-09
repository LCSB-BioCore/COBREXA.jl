"""
Reaction struct.

# Fields
````
id :: String
name :: String
metabolites :: Dict{Metabolite, Float64}
lb :: Float64
ub :: Float64
grr :: Vector{Vector{Gene}}
subsystem :: String
notes :: Dict{String, Vector{String}}
annotation :: Dict{String, Union{Vector{String}, String}}
objective_coefficient :: Float64
````
"""
mutable struct Reaction
    id::String
    name::String
    metabolites::Dict{Metabolite,Float64}
    lb::Float64
    ub::Float64
    grr::Vector{Vector{Gene}}
    subsystem::String
    notes::Dict{String,Vector{String}}
    annotation::Dict{String,Union{Vector{String},String}} # everything is a String[] except sbo, which is a String
    objective_coefficient::Float64

    Reaction(
        id = "";
        name = "",
        metabolites = Dict{Metabolite,Float64}(),
        lb = -1000.0,
        ub = 1000.0,
        grr = Vector{Vector{Gene}}(),
        subsystem = "",
        notes = Dict{String,Vector{String}}(),
        annotation = Dict{String,Union{Vector{String},String}}(),
        objective_coefficient = 0.0,
    ) = new(
        id,
        name,
        metabolites,
        lb,
        ub,
        grr,
        subsystem,
        notes,
        annotation,
        objective_coefficient,
    )

    function Reaction(
        id::String,
        metabolites::Dict{Metabolite,Float64},
        dir = "bidir";
        default_rate = 1000.0,
    )
        if dir == "for"
            lb = 0.0
            ub = default_rate
        elseif dir == "rev"
            lb = -default_rate
            ub = 0.0
        else
            lb = -default_rate
            ub = default_rate
        end
        Reaction(id; metabolites = metabolites, lb = lb, ub = ub)
    end

end
