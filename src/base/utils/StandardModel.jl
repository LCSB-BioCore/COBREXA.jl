
"""
    Base.copy(m::StandardModel)

Shallow copy of a [`StandardModel`](@ref)
"""
Base.copy(m::StandardModel) = StandardModel(
    m.id,
    reactions = m.reactions,
    metabolites = m.metabolites,
    genes = m.genes,
)

"""
    Base.copy(r::Reaction)

Shallow copy of a [`Reaction`](@ref)
"""
Base.copy(r::Reaction) = Reaction(
    r.id;
    metabolites = r.metabolites,
    lb = r.lb,
    ub = r.ub,
    grr = r.grr,
    subsystem = r.subsystem,
    notes = r.notes,
    annotations = r.annotations,
    objective_coefficient = r.objective_coefficient,
)

"""
    Base.copy(m::Metabolite)

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
    Base.copy(g::Gene)

Shallow copy of a [`Gene`](@ref)
"""
Base.copy(g::Gene) = Gene(g.id; notes = g.notes, annotations = g.annotations)

"""
    atom_exchange(model::StandardModel, rxn_id::String)

Specialized version of [`atom_exchange`](@ref) for [`StandardModel`](@ref) that
quickly returns a dictionary mapping the flux-relative stoichiometry of atoms
through a single reaction in `model`.
"""
function atom_exchange(model::StandardModel, rxn_id::String)
    atom_flux = Dict{String,Float64}()
    for (met, stoich_rxn) in model.reactions[rxn_id].metabolites
        adict = metabolite_formula(model, met)
        isnothing(adict) && continue
        for (atom, stoich_molecule) in adict
            atom_flux[atom] = get(atom_flux, atom, 0.0) + stoich_rxn * stoich_molecule
        end
    end
    return atom_flux
end
