"""
buildrxnstring(rxn)

Get rxn in string format for Equilibrator.
"""
function buildrxnstring(rxn::Reaction)
    pos_s = []
    neg_s = []
    for (met, coeff) in rxn.metabolites 
        if coeff > 0.0
            push!(pos_s, "$(coeff) bigg.metabolite:$(met.id[1:end-2])")
        else
            push!(neg_s, "$(abs(coeff)) bigg.metabolite:$(met.id[1:end-2])")    
        end
    end
    return join(neg_s, " + ")*" = "*join(pos_s, " + ") # keep order for ease of use later
end

"""
gibbs_arr = mapGibbs(rxns; dgtype="zero", ph=7.0, ionic_str="100 mM")

Return an dict of rxn.id => ΔG of the specidied dgtype.
"""
function mapGibbs(rxns::Array{Reaction, 1}; dgtype="zero", ph=7.0, ionic_str="100 mM") 
    rxns_strings = [buildrxnstring(rxn) for rxn in rxns]
    if dgtype == "phys"
        py"pygetdg0"(formula, ph, ionic_strength)
        bals, gs, errs = py"pygetdgprimephys"(rxns_strings, ph, ionic_str)
    elseif dgtype == "prime" 
        bals, gs, errs = py"pygetdgprime"(rxns_strings, ph, ionic_str)
    else # "zero"
        bals, gs, errs = py"pygetdg0"(rxns_strings, ph, ionic_str)
    end

    gibbs = Dict{String, Measurement{Float64}}()
    for (rxn, g, err) in zip(rxns, gs, errs)
        if err < 1000.0 # ignore crazy errors
            gibbs[rxn.id] = g ± err 
        else
            gibbs[rxn.id] = 0.0 ± 0.0    
        end
    end
    return gibbs
end

"""
blackbox(fluxres, gibbs, biomassrxn, biomassdg)

Calculate the Gibbs free energy change taking only the external fluxes into account.
"""
function blackbox(fluxres, gibbs, biomassrxn, biomassdg)
    total_ΔG = 0.0 ± 0.0
    for rxnflux in fluxres.rxnfluxes
        if startswith(rxnflux.rxn.id, "EX_")
            if abs(rxnflux.flux) > 1
                total_ΔG -= rxnflux.flux * gibbs[rxnflux.rxn.id] # negative here because formation is not MET -> as used here, but the -> MET 
            end
        end
    end
    total_ΔG += -fluxres.rxnfluxes[fluxres[biomassrxn]].flux*biomassdg # the formation is given correctly here though
    return total_ΔG # units J/gDW/h
end
