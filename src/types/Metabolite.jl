"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)

Default constructor requires the `id` keyword to be assigned.
"""
Base.@kwdef mutable struct Metabolite
    id::String
    name::Maybe{String} = nothing
    formula::Maybe{String} = nothing
    charge::Maybe{Int} = nothing
    compartment::Maybe{String} = nothing
    notes::Notes = Notes()
    annotations::Annotations = Annotations()
end

"""
$(TYPEDSIGNATURES)

A convenience constructor for [`Metabolite`](@ref) taking as first argument the id
of the metabolite. All other kwargs are forwarded to the type constructor.
"""
Metabolite(id; kwargs...) = Metabolite(;id, kwargs...)