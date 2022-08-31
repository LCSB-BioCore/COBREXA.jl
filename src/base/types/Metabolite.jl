"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
mutable struct Metabolite
    id::String
    name::Maybe{String}
    formula::Maybe{String}
    charge::Maybe{Int}
    compartment::Maybe{String}
    notes::Notes
    annotations::Annotations
end

"""
$(TYPEDSIGNATURES)

A constructor for `Metabolite`s.
"""
Metabolite(
    id="";
    name = nothing,
    formula = nothing,
    charge = nothing,
    compartment = nothing,
    notes = Notes(),
    annotations = Annotations(),
) = Metabolite(String(id), name, formula, charge, compartment, notes, annotations)
