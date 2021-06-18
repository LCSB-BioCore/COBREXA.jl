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
    annotations::Annotations
    objective_coefficient::Float64
end

"""
Reaction(
    id = "";
    name = nothing,
    metabolites = Dict{String,Float64}(),
    lb = -_constants.default_reaction_bound,
    ub = _constants.default_reaction_bound,
    grr = nothing,
    subsystem = nothing,
    notes = Notes(),
    annotations = Annotations(),
    objective_coefficient = 0.0,
)

A constructor for Reaction that only takes a reaction `id` and 
assigns default/uninformative values to all the fields that are not
explicitely assigned.
"""
function Reaction(
    id = "";
    name = nothing,
    metabolites = Dict{String,Float64}(),
    lb = -_constants.default_reaction_bound,
    ub = _constants.default_reaction_bound,
    grr = nothing,
    subsystem = nothing,
    notes = Notes(),
    annotations = Annotations(),
    objective_coefficient = 0.0,
)
    mets = Dict(k => float(v) for (k, v) in metabolites)
    return Reaction(
        id,
        name,
        mets,
        lb,
        ub,
        grr,
        subsystem,
        notes,
        annotations,
        objective_coefficient,
    )
end

"""
Reaction(
    id::String,
    metabolites::Dict{String,Union{Int, Float64}},
    dir = :bidirectional;
    default_bound = _constants.default_reaction_bound,
)

Convenience constructor for `Reaction`. The reaction equation is specified using
`metabolites`, which is a dictionary mapping metabolite ids to stoichiometric
coefficients. The direcion of the reaction is set through `dir` which can take
`:bidirectional`, `:forward`, and `:reverse` as values. Finally, the
`default_bound` is the value taken to mean infinity in the context of constraint
based models, often this is set to a very high flux value like 1000.
"""
function Reaction(
    id::String,
    metabolites,
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