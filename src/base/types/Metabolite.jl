"""
$(TYPEDEF)

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

    Metabolite(
        id = "";
        name = nothing,
        formula = nothing,
        charge = nothing,
        compartment = nothing,
        notes = Notes(),
        annotations = Annotations(),
    ) = new(id, name, formula, charge, compartment, notes, annotations)
end
