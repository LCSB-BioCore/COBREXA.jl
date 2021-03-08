"""
    build_rxn_string(rxn::Reaction, compoundtype="kegg")

Get rxn in string format for Equilibrator.
"""
function build_rxn_string(rxn::Reaction, compoundtype="kegg")
    pos_s = []
    neg_s = []

    if compoundtype == "kegg"
        cid = "kegg.compound"
    else
        cid = "bigg.metabolite"
    end

    for (met, coeff) in rxn.metabolites
        metid = get(met.annotation, cid, [""])[1]
        metid == "" && continue 
        if coeff > 0.0
            if compoundtype == "kegg" 
                push!(pos_s, "$(coeff) KEGG:$metid")
            else
                push!(pos_s, "$(coeff) bigg.metabolite:$metid")
            end
        else
            if compoundtype == "kegg" 
                push!(neg_s, "$(abs(coeff)) KEGG:$metid")
            else
                push!(neg_s, "$(abs(coeff)) bigg.metabolite:$metid")
            end 
        end
    end
    return join(neg_s, " + ")*" = "*join(pos_s, " + ") # keep order for ease of use later
end

"""
    map_gibbs_rxns(rxns::Array{Reaction, 1}; dgtype="zero", ph=7.0, ionic_str="100 mM", usekegg=true)

Return a dict of rxn.id => ΔG of the specified dgtype.
Takes as inputs an array of reactions, `rxns`, and optional keyword arguments `dgtype`, which specifies which type of ΔG calculation should be returned.
Valid options for `dgtype` are "phys" and "prime" (anything else is "zero"). 
Ionic strength can be set through `ionic_str` which takes a string input, e.g. "150 mM".
Only BIGG and KEGG metabolite identifiers are supported, i.e. the reaction needs to have a KEGG or BIGG `id` listed in the annotation field in the reaction struct.
By default KEGG annotations are used to build the reaction strings that are fed to Equilibrator. Note that the first metabolite `id` is used.
"""
function map_gibbs_rxns(rxns::Array{Reaction, 1}; dgtype="zero", ph=7.0, ionic_str="100 mM", usekegg=true) 
    if usekegg
        rxns_strings = [build_rxn_string(rxn, "kegg") for rxn in rxns]
    else
        rxns_strings = [build_rxn_string(rxn, "bigg") for rxn in rxns]
    end
    
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

    balances = Dict{String, Float64}()
    for (rxn, bal) in zip(rxns, bals)
        balances[rxn.id] = bal
    end

    return gibbs, balances
end

"""
    map_gibbs_external(fluxres::Dict{String, Float64}, gibbs)

Calculate the Gibbs free energy change taking only the external fluxes into account.
NB: you need to account for the biomass function separately.

Fluxres can be both a ReactionFluxes object or a Dict with rxnid -> flux.
"""
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
    map_gibbs_internal(fluxres, gibbs, biomassid="BIOMASS")

Calculate the Gibbs free energy change taking only the internal fluxes into account.
NB: you need to account for the biomass function separately. NB: 

Fluxres are a Dict with rxnid -> flux.
"""
function map_gibbs_internal(fluxres::Dict{String, Float64}, gibbs, biomassid="BIOMASS")
    total_ΔG = 0.0 ± 0.0
    missing_flux = 0.0
    found_flux = 0.0
    for (rxnid, v) in fluxres
        if !startswith(rxnid, "EX_") && !contains(rxnid, biomassid) # ignore exchange reactions and biomass function
            if gibbs[rxnid] ≈ 0.0
                missing_flux += abs(v)
            else
                found_flux += abs(v)
            end 
            total_ΔG += v * gibbs[rxnid] # add because this is not formation but rather just adding equations (the flux direction sign compensates)
        end
    end
    return total_ΔG, missing_flux/(missing_flux+found_flux) # units J/gDW/h
end