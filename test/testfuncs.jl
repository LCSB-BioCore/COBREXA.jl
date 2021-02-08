"""
Shouldn't test the identity of S, b, lb, ub because the storage format may differ i.e. order of rxns
and mets could be different. Rather test indirectly.
"""
function model_comparison_test(model1, model2)
    # test if blank model is given - automatic fail
    (isempty(model1.rxns) || isempty(model2.rxns)) ? (return false) : nothing

    # test same rxn and met ids
    rxns1 = [r.id for r in model1.rxns]
    rxns2 = [r.id for r in model2.rxns]
    mets1 = [m.id for m in model1.mets]
    mets2 = [m.id for m in model2.mets]
    rxn_diff = (isempty(setdiff(rxns1, rxns2)) && isempty(setdiff(rxns2, rxns1))) ? true : false
    met_diff = (isempty(setdiff(mets1, mets2)) && isempty(setdiff(mets2, mets1))) ? true : false

    # test lb and ub the same (indirectly)
    lbs1 = [r.lb for r in model1.rxns]
    lbs2 = [r.lb for r in model2.rxns]
    ubs1 = [r.lb for r in model1.rxns]
    ubs2 = [r.lb for r in model2.rxns]

    lb_same = sum(abs, lbs1) == sum(abs, lbs2)
    ub_same = sum(abs, ubs1) == sum(abs, ubs2)
        
    # test grrs
    model1_grr_keys = keys(model1.grrs)
    model2_grr_keys = keys(model2.grrs)
    same_grr_keys = (isempty(setdiff(model1_grr_keys, model2_grr_keys)) && isempty(setdiff(model2_grr_keys, model1_grr_keys))) ? true : false
    same_grr = true # change if inconsistency detected
    if same_grr_keys
        grr_keys = keys(model1.grrs)
        for grr_key in grr_keys
            model1_gene_sets = model1.grrs[grr_key]
            model2_gene_sets = model2.grrs[grr_key]
            for model1_gene_set in model1_gene_sets
                found = false
                for model2_gene_set in model2_gene_sets
                    found = (isempty(setdiff(model1_gene_set, model2_gene_set)) && isempty(setdiff(model2_gene_set, model1_gene_set)))
                    found && break
                end
                !found && (same_grr = false; break)
            end
            !same_grr && break
        end
    else
        same_grr = false
    end

    return all([rxn_diff, met_diff, lb_same, ub_same, same_grr])
end

function read_write_read_test(model, format)
    tmpfile = "temp."*format

    CobraTools.savemodel(model, tmpfile)
    tmpmodel = CobraTools.readmodel(tmpfile)

    rm(tmpfile)
    model_comparison_test(model, tmpmodel)
end

function rxn_construction_test(model)
    # note, this test is for iJO1366
    rxn_original = findfirst(model.rxns, "NADH16pp")
    nadh = findfirst(model.mets, "nadh_c")
    h_c = findfirst(model.mets, "h_c")
    q8 = findfirst(model.mets, "q8_c")
    q8h2 = findfirst(model.mets, "q8h2_c")
    nad = findfirst(model.mets, "nad_c")
    h_p = findfirst(model.mets, "h_p")
    
    rxn = 1.0*nadh + 4.0*h_c + 1.0*q8 ⟶  1.0*q8h2 + 1.0*nad + 3.0*h_p
    check_bounds_forward = (rxn.lb == 0.0 && rxn.ub > 0.0) ? true : false

    rxn = 1.0*nadh + 4.0*h_c + 1.0*q8 ← 1.0*q8h2 + 1.0*nad + 3.0*h_p
    check_bounds_reverse = (rxn.lb < 0.0 && rxn.ub == 0.0) ? true : false
    
    rxn = 1.0*nadh + 4.0*h_c + 1.0*q8 ↔ 1.0*q8h2 + 1.0*nad + 3.0*h_p
    check_bounds_bidir = (rxn.lb < 0.0 && rxn.ub > 0.0) ? true : false
    
    rxn = 1.0*nadh → ∅
    check_ex_out = length(rxn.metabolites) == 1 ? true : false 
    rxn = ∅ → nadh
    check_ex_in = length(rxn.metabolites) == 1 ? true : false 
    
    rxn = 1.0*nadh + 4.0*h_c + 1.0*q8 ⟶  1.0*q8h2 + 1.0*nad + 3.0*h_p
    rxn_mets_coeffs = prod(values(rxn.metabolites)) == -12 ? true : false
    rxn_mets = ("q8h2_c" in [x.id for x  in keys(rxn.metabolites)]) ? true : false # getting one right suggests it works

    return all([check_bounds_forward, check_bounds_reverse, check_bounds_bidir, check_ex_out, check_ex_in, rxn_mets, rxn_mets_coeffs])
end

function fba_test(model)    
    biomass_rxn = findfirst(model.rxns, "BIOMASS_Ec_iJO1366_WT_53p95M")
    solobj = CobraTools.fba(model, biomass_rxn)
    return solobj.obj ≈ 0.9865144469529787
end

function pfba_test(model)    
    biomass_rxn = findfirst(model.rxns, "BIOMASS_Ec_iJO1366_WT_53p95M")
    solobj = CobraTools.pfba(model, biomass_rxn)
    return solobj.obj ≈ 15546.145490407944
end