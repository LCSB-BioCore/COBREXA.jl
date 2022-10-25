"""
$(TYPEDSIGNATURES)

Return two dictionaries of metabolite `id`s mapped to reactions that consume or
produce them, given the flux distribution supplied in `flux_dict`.
"""
function metabolite_fluxes(model::AbstractMetabolicModel, flux_dict::Dict{String,Float64})
    S = stoichiometry(model)
    rids = reactions(model)
    mids = metabolites(model)

    producing = Dict{String,Dict{String,Float64}}()
    consuming = Dict{String,Dict{String,Float64}}()
    for (row, mid) in enumerate(mids)
        for (col, rid) in enumerate(rids)
            mf = flux_dict[rid] * S[row, col]
            if mf < 0 # consuming rxn
                if haskey(consuming, mid)
                    consuming[mid][rid] = mf
                else
                    consuming[mid] = Dict(rid => mf)
                end
            elseif mf >= 0 # producing rxn
                if haskey(producing, mid)
                    producing[mid][rid] = mf
                else
                    producing[mid] = Dict(rid => mf)
                end
            end
        end
    end
    return consuming, producing
end

"""
$(TYPEDSIGNATURES)

Return a dictionary mapping the flux of atoms across a flux solution given by
`reaction_fluxes` using the reactions in `model` to determine the appropriate stoichiometry.

Note, this function ignores metabolites with no formula assigned to them, no error message
will be displayed.

Note, if a model is mass balanced there should be not net flux of any atom. By removing
reactions from the flux_solution you can investigate how that impacts the mass balances.

# Example
```
# Find flux of Carbon through all metabolic reactions except the biomass reaction
delete!(fluxes, "BIOMASS_Ecoli_core_w_GAM")
atom_fluxes(model, fluxes)["C"]
```
"""
function atom_fluxes(model::AbstractMetabolicModel, reaction_fluxes::Dict{String,Float64})
    rids = reactions(model)
    atom_flux = Dict{String,Float64}()
    for (ridx, rid) in enumerate(rids)
        haskey(reaction_fluxes, rid) || continue
        rflux = reaction_fluxes[rid]
        for (mid, mstoi) in reaction_stoichiometry(model, rid)
            atoms = metabolite_formula(model, mid)
            isnothing(atoms) && continue # missing formulas are ignored
            for (atom, abundance) in atoms
                atom_flux[atom] = get(atom_flux, atom, 0.0) + rflux * mstoi * abundance
            end
        end
    end
    return atom_flux
end
