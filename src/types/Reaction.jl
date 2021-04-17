"""
Reaction struct.

Note that the `metabolite` field maps metabolite `id`s to stoichiometric coefficients.

# Fields
````
id :: String
name :: String
metabolites :: Dict{String, Float64}
lb :: Float64
ub :: Float64
grr :: Vector{Vector{String}}
subsystem :: String
notes :: Dict{String, Vector{String}}
annotation :: Dict{String, Union{Vector{String}, String}}
objective_coefficient :: Float64
````
"""
mutable struct Reaction
    id::String
    name::String
    metabolites::Dict{String,Float64} # reaction id => stoichiometric coefficient
    lb::Float64
    ub::Float64
    grr::Vector{Vector{String}} # [[rxn_id1, rxn_id2], [rxn_id3]]
    subsystem::String
    notes::Dict{String,Vector{String}}
    annotation::Dict{String,Union{Vector{String},String}} # everything is a String[] except sbo, which is a String
    objective_coefficient::Float64

    Reaction(
        id = "";
        name = "",
        metabolites = Dict{String,Float64}(),
        lb = -_constants.default_reaction_bound,
        ub = _constants.default_reaction_bound,
        grr = Vector{Vector{String}}(),
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
