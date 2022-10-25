"""
$(TYPEDEF)

A structure for representing a single reaction in a [`ObjectModel`](@ref).

# Fields
$(TYPEDFIELDS)
"""
mutable struct Reaction
    id::String
    name::Maybe{String}
    metabolites::Dict{String,Float64}
    lower_bound::Float64
    upper_bound::Float64
    grr::Maybe{GeneAssociation}
    subsystem::Maybe{String}
    notes::Notes
    annotations::Annotations
    objective_coefficient::Float64
end

"""
$(TYPEDSIGNATURES)

A constructor for Reaction that only takes a reaction `id` and
assigns default/uninformative values to all the fields that are not
explicitely assigned.
"""
function Reaction(
    id = "";
    name = nothing,
    metabolites = Dict{String,Float64}(),
    lower_bound = -constants.default_reaction_bound,
    upper_bound = constants.default_reaction_bound,
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
        lower_bound,
        upper_bound,
        grr,
        subsystem,
        notes,
        annotations,
        objective_coefficient,
    )
end

"""
$(TYPEDSIGNATURES)

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
    default_bound = constants.default_reaction_bound,
)
    if dir == :forward
        lower_bound = 0.0
        upper_bound = default_bound
    elseif dir == :reverse
        lower_bound = -default_bound
        upper_bound = 0.0
    elseif dir == :bidirectional
        lower_bound = -default_bound
        upper_bound = default_bound
    else
        throw(DomainError(dir, "unsupported direction"))
    end
    Reaction(
        id;
        metabolites = metabolites,
        lower_bound = lower_bound,
        upper_bound = upper_bound,
    )
end
