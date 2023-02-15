
"""
$(TYPEDSIGNATURES)

Shallow copy of a [`ObjectModel`](@ref)
"""
Base.copy(m::ObjectModel) = ObjectModel(
    reactions = m.reactions,
    metabolites = m.metabolites,
    genes = m.genes,
    objective = m.objective,
    notes = m.notes,
    annotations = m.annotations,
)

"""
$(TYPEDSIGNATURES)

Shallow copy of a [`Reaction`](@ref)
"""
Base.copy(r::Reaction) = Reaction(
    r.id;
    metabolites = r.metabolites,
    lower_bound = r.lower_bound,
    upper_bound = r.upper_bound,
    gene_associations = r.gene_associations,
    subsystem = r.subsystem,
    notes = r.notes,
    annotations = r.annotations,
)

"""
$(TYPEDSIGNATURES)

Shallow copy of a [`Metabolite`](@ref)
"""
Base.copy(m::Metabolite) = Metabolite(
    m.id;
    formula = m.formula,
    charge = m.charge,
    compartment = m.compartment,
    notes = m.notes,
    annotations = m.annotations,
)

"""
$(TYPEDSIGNATURES)

Shallow copy of a [`Gene`](@ref)
"""
Base.copy(g::Gene) = Gene(g.id; notes = g.notes, annotations = g.annotations)

"""
$(TYPEDSIGNATURES)

Return the lower bounds for all reactions in `model`.
Order matches that of the reaction IDs returned by [`reactions`](@ref).
"""
lower_bounds(model::ObjectModel)::Vector{Float64} =
    [model.reactions[rxn].lower_bound for rxn in reactions(model)]

"""
$(TYPEDSIGNATURES)

Return the upper bounds for all reactions in `model`.
Order matches that of the reaction IDs returned in [`reactions`](@ref).
"""
upper_bounds(model::ObjectModel)::Vector{Float64} =
    [model.reactions[rxn].upper_bound for rxn in reactions(model)]
