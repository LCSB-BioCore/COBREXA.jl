"""
Shouldn't test the identity of S, b, lb, ub because the storage format may differ i.e. order of rxns
and mets could be different. Rather test indirectly.
"""
function model_comparison_test(model1, model2)
    # test if blank model is given - automatic fail
    (isempty(model1.reactions) || isempty(model2.reactions)) ? (return false) : nothing

    # test same rxn and met ids
    rxns1 = [r.id for r in model1.reactions]
    rxns2 = [r.id for r in model2.reactions]
    mets1 = [m.id for m in model1.metabolites]
    mets2 = [m.id for m in model2.metabolites]
    rxns_same = (isempty(setdiff(rxns1, rxns2)) && isempty(setdiff(rxns2, rxns1))) ? true : false
    mets_same = (isempty(setdiff(mets1, mets2)) && isempty(setdiff(mets2, mets1))) ? true : false

    # test lb and ub the same (indirectly)
    lbs1 = [r.lb for r in model1.reactions]
    lbs2 = [r.lb for r in model2.reactions]
    ubs1 = [r.lb for r in model1.reactions]
    ubs2 = [r.lb for r in model2.reactions]

    lbs_same = sum(abs, lbs1) == sum(abs, lbs2)
    ubs_same = sum(abs, ubs1) == sum(abs, ubs2)
        
    S1, _, _, _ = get_core_model(model1)
    S2, _, _, _ = get_core_model(model2)
    Ss_same = sum(S1) == sum(S2) ? true : false

    grrs_same = true
    for rxn1 in model1.reactions
        rxn2 = findfirst(model2.reactions, rxn1.id)
        rxn1_grr_string = string(split(CobraTools.unparse_grr(rxn1.grr), "")) # use hash
        rxn2_grr_string = string(split(CobraTools.unparse_grr(rxn2.grr), ""))
        rxn1_ints = sum(Int.([x[1] for x in rxn1_grr_string]))
        rxn2_ints = sum(Int.([x[1] for x in rxn2_grr_string]))
        if rxn1_ints != rxn2_ints
            grrs_same = false
            break
        end
    end

    return all([rxns_same, mets_same, lbs_same, ubs_same, Ss_same, grrs_same])
end

"""
Test if model is the same after it was read in, saved, and then re-read. 
"""
function read_write_read_test(model, format)
    tmpfile = joinpath("data", "temp."*format)

    CobraTools.save_model(model, tmpfile)
    tmpmodel = CobraTools.read_model(tmpfile)

    rm(tmpfile)
    model_comparison_test(model, tmpmodel)
end



# function atom_test(model)
#     biomass_rxn = findfirst(model.rxns, "BIOMASS_Ec_iJO1366_WT_53p95M")
#     pfbasol = CobraTools.pfba(model, biomass_rxn)
#     ad = CobraTools.atom_exchange(pfbasol)
#     return ad["C"]/ad["H"] ≈ 0.6362486422376349
# end