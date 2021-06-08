
"""
    Base.copy(m::StandardModel)

Shallow copy of a [`StandardModel`](@ref)
"""
Base.copy(m::StandardModel) = StandardModel(m.id, m.reactions, m.metabolites, m.genes)

"""
    Base.copy(r::Reaction)

Shallow copy of a [`Reaction`](@ref)
"""
Base.copy(r::Reaction) = Reaction(
    r.id;
    name = r.name,
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
    name = m.name,
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
Base.copy(g::Gene) = Gene(g.id; name = g.name, notes = g.notes, annotations = g.annotations)

"""
    atom_exchange(flux_dict::Dict{String, Float64}, model::StandardModel)

Return a dictionary mapping the flux of atoms across the boundary of the model 
given `flux_dict` (the solution of a constraint based analysis) of reactions in `model`.
"""
function atom_exchange(flux_dict::Dict{String,Float64}, model::StandardModel)
    atom_flux = Dict{String,Float64}()
    for (rxn_id, flux) in flux_dict
        if is_boundary(model.reactions[rxn_id])
            for (met, stoich_rxn) in model.reactions[rxn_id].metabolites
                adict = get_atoms(model.metabolites[met])
                for (atom, stoich_molecule) in adict
                    atom_flux[atom] = get(atom_flux, atom, 0.0) + flux * stoich_rxn * stoich_molecule
                end
            end
        end
    end
    return atom_flux
end

"""
    atom_exchange(rxn_id::String, model::StandardModel)

Return a dictionary mapping the flux of atoms through a reaction in `model`.
"""
function atom_exchange(rxn_id::String, model::StandardModel)
    atom_flux = Dict{String,Float64}()
    for (met, stoich_rxn) in model.reactions[rxn_id].metabolites
        adict = get_atoms(model.metabolites[met])
        for (atom, stoich_molecule) in adict
            atom_flux[atom] = get(atom_flux, atom, 0.0) + stoich_rxn * stoich_molecule
        end
    end
    return atom_flux
end

"""
    metabolite_fluxes(flux_dict::Dict{String, Float64}, model::StandardModel)

Return two dictionaries of metabolite `id`s mapped to reactions that consume or 
produce them given the flux distribution supplied in `fluxdict`.
"""
function metabolite_fluxes(flux_dict::Dict{String,Float64}, model::StandardModel)
    S = stoichiometry(model)
    rxnids = reactions(model)
    metids = metabolites(model)

    producing = Dict{String,Dict{String,Float64}}()
    consuming = Dict{String,Dict{String,Float64}}()
    for (row, metid) in enumerate(metids)
        for (col, rxnid) in enumerate(rxnids)
            mf = flux_dict[rxnid] * S[row, col]
            # ignore zero flux
            if mf < -_constants.tolerance # consuming rxn
                if haskey(consuming, metid)
                    consuming[metid][rxnid] = mf
                else
                    consuming[metid] = Dict(rxnid => mf)
                end
            elseif mf > _constants.tolerance
                if haskey(producing, metid)
                    producing[metid][rxnid] = mf
                else
                    producing[metid] = Dict(rxnid => mf)
                end
            end
        end
    end
    return consuming, producing
end
