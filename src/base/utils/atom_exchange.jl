"""
    atom_exchange(model::MetabolicModel, flux_dict::Dlicict{String, Float64})

Return a dictionary mapping the flux of atoms across the boundary of the model 
given `flux_dict` (the solution of a constraint based analysis) of reactions in `model`.
"""
function atom_exchange(model::MetabolicModel, flux_dict::Dict{String,Float64})
    S = stoichiometry(model)
    mids = metabolites(model)
    rids = reactions(model)

    atom_flux = Dict{String,Float64}()
    for (ridx, rid) in enumerate(rids)
        rflux = flux_dict[rid]
        for (midx, mstoi) in zip(findnz(S[:, ridx])...)
            atoms = metabolite_formula(model, mids[midx])
            isnothing(atoms) && continue
            for (atom, abundance) in atoms
                atom_flux[atom] = get(atom_flux, atom, 0.0) + rflux * mstoi * abundance
            end
        end
    end
    return atom_flux
end
