"""
Shouldn't test the identity of S, b, lb, ub because the storage format may differ i.e. order of rxns
and mets could be different. Rather test indirectly.
"""
function model_comparison(model1, model2)
    # test if blank model is given - automatic fail
    (isempty(model1.coremodel.S) || isempty(model2.coremodel.S)) ? (return false) : nothing

    # test same rxn and met ids
    rxns1 = [r.id for r in model1.rxns]
    rxns2 = [r.id for r in model2.rxns]
    mets1 = [m.id for m in model1.mets]
    mets2 = [m.id for m in model2.mets]
    rxn_diff = (isempty(setdiff(rxns1, rxns2)) && isempty(setdiff(rxns2, rxns1))) ? true : false
    met_diff = (isempty(setdiff(mets1, mets2)) && isempty(setdiff(mets2, mets1))) ? true : false

    # test same S and b shapes
    S_size = all(size(model1.coremodel.S) .== size(model2.coremodel.S))
    b_size = all(size(model1.coremodel.b) .== size(model2.coremodel.b))

    # test lb and ub the same (indirectly)
    lbs1 = [r.lb for r in model1.rxns]
    lbs2 = [r.lb for r in model2.rxns]
    ubs1 = [r.lb for r in model1.rxns]
    ubs2 = [r.lb for r in model2.rxns]

    lb_same = sum(abs, lbs1) == sum(abs, lbs2)
    ub_same = sum(abs, ubs1) == sum(abs, ubs2)
    
    # test if S and b the same (indirectly)
    if S_size
        input_vector = rand(size(model1.coremodel.S, 2))
        S_same = sum(abs, model1.coremodel.S * input_vector) ≈ sum(abs, model2.coremodel.S * input_vector)
    else
        S_same = false
    end
    b_same = sum(abs, model1.coremodel.b) == sum(abs, model2.coremodel.b)
    
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

    return all([rxn_diff, met_diff, S_size, b_size, lb_same, ub_same, S_same, b_same, same_grr])
end

function read_write_read(model, format)
    tmpfile = "temp."*format

    CobraTools.savemodel(model, tmpfile)
    tmpmodel = CobraTools.readmodel(tmpfile)

    rm(tmpfile)
    model_comparison(model, tmpmodel)
end