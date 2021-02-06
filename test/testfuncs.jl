"""
Shouldn't test the identity of S, b, lb, ub because the storage format may differ i.e. order of rxns
and mets could be different. Rather test indirectly.
"""
function model_comparison(model1, model2)
    # test if blank model is given - automatic fail
    (isempty(model1.S) || isempty(model2.S)) ? (return false) : nothing

    # test same rxn and met ids
    rxn_diff = (isempty(setdiff(model1.rxns, model2.rxns)) && isempty(setdiff(model2.rxns, model1.rxns))) ? true : false
    met_diff = (isempty(setdiff(model1.mets, model2.mets)) && isempty(setdiff(model2.mets, model1.mets))) ? true : false

    # test same S and b shapes
    S_size = all(size(model1.S) .== size(model2.S))
    b_size = all(size(model1.b) .== size(model2.b))

    # test lb and ub the same (indirectly)
    lb_same = sum(abs, model1.lbs) == sum(abs, model2.lbs)
    ub_same = sum(abs, model1.ubs) == sum(abs, model2.ubs)
    
    # test if S and b the same (indirectly)
    if S_size
        input_vector = rand(size(model1.S, 2))
        S_same = sum(abs, model1.S * input_vector) ≈ sum(abs, model2.S * input_vector)
    else
        S_same = false
    end
    b_same = sum(abs, model1.b) == sum(abs, model2.b)
    
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

