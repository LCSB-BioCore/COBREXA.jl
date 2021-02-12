"""
buildrxnstring(rxn)

Get rxn in string format for Equilibrator.
"""
function build_rxn_string(rxn::Reaction)
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
function map_gibbs_rxns(rxns::Array{Reaction, 1}; dgtype="zero", ph=7.0, ionic_str="100 mM") 
    rxns_strings = [build_rxn_string(rxn) for rxn in rxns]
    if dgtype == "phys"
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
map_gibbs_external(fluxres, gibbs)

Calculate the Gibbs free energy change taking only the external fluxes into account.
NB: you need to account for the biomass function separately.

Fluxres can be both a ReactionFluxes object or a Dict with rxnid -> flux.
"""
function map_gibbs_external(fluxres::ReactionFluxes, gibbs)
    total_ΔG = 0.0 ± 0.0
    missing_flux = 0.0
    for (i, rxn) in enumerate(fluxres.rxns)
        if startswith(rxn.id, "EX_")
            if gibbs[rxn.id] ≈ 0.0
                missing_flux += abs(fluxres.fluxes[i])
            end 
            total_ΔG -= fluxres.fluxes[i] * gibbs[rxn.id] # negative here because formation is not MET -> as used here, but the -> MET 
        end
    end
    return total_ΔG, missing_flux/sum(abs, fluxres.fluxes) # units J/gDW/h
end

function map_gibbs_external(fluxres::Dict{String, Float64}, gibbs)
    total_ΔG = 0.0 ± 0.0
    missing_flux = 0.0
    for (rxnid, v) in fluxres
        if startswith(rxnid, "EX_")
            if gibbs[rxnid] ≈ 0.0
                missing_flux += abs(v)
            end    
            total_ΔG -= v * gibbs[rxnid] # negative here because "combustion" is actually Gibbs value not formation  
        end
    end
    return total_ΔG, missing_flux/sum(abs, values(fluxres)) # units J/gDW/h
end

"""
map_gibbs_internal(fluxres, gibbs)

Calculate the Gibbs free energy change taking only the internal fluxes into account.
NB: you need to account for the biomass function separately

Fluxres can be both a ReactionFluxes object or a Dict with rxnid -> flux.
"""
function map_gibbs_internal(fluxres::ReactionFluxes, gibbs)
    total_ΔG = 0.0 ± 0.0
    missing_flux = 0.0
    for (i, rxn) in enumerate(fluxres.rxns)
        if !startswith(rxn.id, "EX_") # ignore exchange reactions
            if abs(fluxres.fluxes[i]) > 1e-8
                if gibbs[rxn.id] ≈ 0.0
                    missing_flux += abs(fluxres.fluxes[i])
                end
                total_ΔG += fluxres.fluxes[i] * gibbs[rxn.id] # add because this is not formation but rather just adding equations (the flux direction sign compensates) 
            end
        end
    end
    return total_ΔG, missing_flux/sum(abs, fluxres.fluxes) # units J/gDW/h
end

function map_gibbs_internal(fluxres::Dict{String, Float64}, gibbs)
    total_ΔG = 0.0 ± 0.0
    missing_flux = 0.0
    for (rxnid, v) in fluxres
        if !startswith(rxnid, "EX_") # ignore exchange reactions 
            if gibbs[rxnid] ≈ 0.0
                missing_flux += abs(v)
            end 
            total_ΔG += v * gibbs[rxnid] # add because this is not formation but rather just adding equations (the flux direction sign compensates)
        end
    end
    return total_ΔG, missing_flux/sum(abs, values(fluxres)) # units J/gDW/h
end

"""
save_Gibbs(path, gibbs)

Save Gibbs dict. Saved as String => [mag, err]
"""
function save_Gibbs(path, gibbs)
    decomp = Dict{String, Array{Float64,1}}()
    for (k, v) in gibbs
        decomp[k] = [Measurements.value(v), Measurements.uncertainty(v)]
    end
    JLD.save(path, "gibbs", decomp)
end

"""
load_Gibbs(path)

Load Gibbs dict. Loads String => [mag, err]
"""
function load_Gibbs(path)
    decomp = JLD.load(path, "gibbs")
    gibbs = Dict{String, Measurement{Float64}}()
    for (k, vs) in decomp
        gibbs[k] = vs[1] ± vs[2]
    end
    return gibbs
end