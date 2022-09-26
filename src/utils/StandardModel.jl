
"""
$(TYPEDSIGNATURES)

Shallow copy of a [`StandardModel`](@ref)
"""
Base.copy(m::StandardModel) = StandardModel(
    m.id,
    reactions = m.reactions,
    metabolites = m.metabolites,
    genes = m.genes,
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
    grr = r.grr,
    subsystem = r.subsystem,
    notes = r.notes,
    annotations = r.annotations,
    objective_coefficient = r.objective_coefficient,
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
