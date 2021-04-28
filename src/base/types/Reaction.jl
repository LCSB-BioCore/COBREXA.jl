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
    name::Maybe{String}
    metabolites::Dict{String,Float64}
    lb::Float64
    ub::Float64
    grr::Maybe{GeneAssociation}
    subsystem::Maybe{String}
    notes::Notes
    annotation::Annotations# everything is a String[]
    objective_coefficient::Float64 # defaults to 0.0

    Reaction(
        id = "";
        name = nothing,
        metabolites = Dict{String,Float64}(),
        lb = -_constants.default_reaction_bound,
        ub = _constants.default_reaction_bound,
        grr = nothing,
        subsystem = nothing,
        notes = Dict{String,Vector{String}}(),
        annotations = Dict{String,Vector{String}}(),
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
        annotations,
        objective_coefficient,
    )

    # this constructor is only used internally
    function Reaction(
        id::String,
        metabolites::Dict{String,Float64},
        dir = :bidirectional;
        default_bound = _constants.default_reaction_bound,
    )
        if dir == :forward
            lb = 0.0
            ub = default_bound
        elseif dir == :reverse
            lb = -default_bound
            ub = 0.0
        else
            lb = -default_bound
            ub = default_bound
        end
        Reaction(id; metabolites = metabolites, lb = lb, ub = ub)
    end
end
