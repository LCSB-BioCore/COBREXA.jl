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
    tmpfile = "temp."*format

    CobraTools.save_model(model, tmpfile)
    tmpmodel = CobraTools.read_model(tmpfile)

    rm(tmpfile)
    model_comparison_test(model, tmpmodel)
end

"""
Test if gene can be constructed and manipulated.
"""
function test_gene()
    g = Gene()
    g.id = "gene1"
    g.name = "gene_name"
    g.notes = Dict("notes"=>["blah", "blah"])
    g.annotation = Dict("sboterm" => "sbo", "ncbigene" => ["ads", "asds"])

    println(g) # test IO
    
    g2 = Gene("gene2")
    
    genes = [g, g2]

    println(genes) # test IO
    
    ind = genes[g]
    if ind != 1
        return false
    end
    
    gg = findfirst(genes, g2.id)
    if gg.id != g2.id
        return false
    end
    
    g3 = Gene("g3")
    g3.annotation = Dict("ncbigene" => "sbo", "ncbigene" => ["ads", "asds"])
    
    dup, ind = check_duplicate_annotations(genes, g3)
    if !dup || ind != 1 
        return false
    end

    return true
end

"""
Test if metabolite can be constructed and manipulated.
"""
function test_metabolite()
    m1 = Metabolite()
    m1.id = "met1"
    m1.name = "metabolite 1"
    m1.formula = "C6H12O6N"
    m1.charge = 1
    m1.compartment = "c"
    m1.notes = Dict("notes"=>["blah", "blah"])
    m1.annotation = Dict("sboterm" => "sbo", "kegg.compound" => ["ads", "asds"])
    
    println(m1) # test IO

    m2 = Metabolite("met2")
    m2.formula = "C6H12O6N"
    
    m3 = Metabolite("met3")
    m3.formula = "X"
    m3.annotation = Dict("sboterm" => "sbo", "kegg.compound" => ["ad2s", "asds"])
    
    mets = [m1, m2, m3]

    println(mets) # test IO
    
    ind = mets[m2]
    if ind != 2
        return false
    end
    
    mm = findfirst(mets, "met3")
    if mm.id != m3.id
        return false
    end
    
    dup, ind = check_duplicate_annotations(mets, m3)
    if !dup || ind !=3
        return false
    end
    
    mms = check_same_formula([m3, m1], m2)
    if length(mms) != 1
        return false
    end
    
    ats = get_atoms(mms[1])
    if ats["C"] != 6 && ats["N"] != 1
        return false
    end
    
    return true
end

"""
Test if reaction can be constructed and manipulated.
"""
function test_reaction()
    m1 = Metabolite("m1")
    m1.formula = "C2H3"
    m2 = Metabolite("m2")
    m2.formula = "H3C2"
    m3 = Metabolite("m3")
    m4 = Metabolite("m4")
    
    g1 = Gene("g1")
    g2 = Gene("g2")
    g3 = Gene("g3")
    
    r1 = Reaction()
    r1.id = "r1"
    r1.name = "reaction 1"
    r1.metabolites = Dict(m1 => -1.0, m2 => 1.0)
    r1.lb = -100.0
    r1.ub = 100.0
    r1.grr = [[g1, g2], [g3]]
    r1.subsystem = "glycolysis"
    r1.notes = Dict("notes"=>["blah", "blah"])
    r1.annotation = Dict("sboterm" => "sbo", "biocyc" => ["ads", "asds"])
    r1.objective_coefficient = 1.0

    println(r1) # test IO
    
    r2 = Reaction("r2", Dict(m1 => -2.0, m4 => 1.0), "rev")
    if r2.lb != -1000.0 && r2.ub != 0.0
        return false
    end
    
    r3 = Reaction("r3", Dict(m3 => -1.0, m4 => 1.0), "for")
    if r3.lb != 0.0 && r3.ub != 1000.0
        return false
    end
    
    rxns = [r1, r2, r3]

    println(rxns) # test IO
    
    ind = rxns[r3]
    if ind != 3
        return false
    end
    
    rr = findfirst(rxns, "r2")
    if rr.id != r2.id
        return false
    end
    
    r4 = Reaction("r4", Dict(m3 => -1.0, m4 => 1.0), "bidir")
    r4.annotation = Dict("sboterm" => "sbo", "biocyc" => ["ads", "asds"])
    if r4.lb != -1000.0 && r4.ub != 1000.0
        return false
    end
    
    dup, ind = check_duplicate_annotations(rxns, r4)
    if !dup || ind !=1
        return false
    end
    
    dup, ind = check_duplicate_reaction(rxns, r4)
    if !dup || ind !=3
        return false
    end
    
    bal, d = is_mass_balanced(r1)
    if !bal
        return false
    end

    return true
end

"""
Test if model can be constructed and manipulated (basic).
"""
function test_model()
    m1 = Metabolite("m1")
    m1.formula = "C2H3"
    m2 = Metabolite("m2")
    m2.formula = "H3C2"
    m3 = Metabolite("m3")
    m4 = Metabolite("m4")

    g1 = Gene("g1")
    g2 = Gene("g2")
    g3 = Gene("g3")

    r1 = Reaction()
    r1.id = "r1"
    r1.name = "reaction 1"
    r1.metabolites = Dict(m1 => -1.0, m2 => 1.0)
    r1.lb = -100.0
    r1.ub = 100.0
    r1.grr = [[g1, g2], [g3]]
    r1.subsystem = "glycolysis"
    r1.notes = Dict("notes"=>["blah", "blah"])
    r1.annotation = Dict("sboterm" => "sbo", "biocyc" => ["ads", "asds"])
    r1.objective_coefficient = 1.0

    r2 = Reaction("r2", Dict(m1 => -2.0, m4 => 1.0), "rev")
    r3 = Reaction("r3", Dict(m3 => -1.0, m4 => 1.0), "for")
    r4 = Reaction("r4", Dict(m3 => -1.0, m4 => 1.0), "bidir")
    r4.annotation = Dict("sboterm" => "sbo", "biocyc" => ["ads", "asds"])

    mets = [m1, m2, m3, m4]
    genes = [g1, g2, g3]
    rxns = [r1, r2, r3, r4]
    
    model = Model()
    model.id = "model"
    model.reactions = rxns
    model.metabolites = mets
    model.genes = genes
    
    println(model) # test IO

    ind = model[r2]
    if ind != 2
        return false
    end

    ind = model[g3]
    if ind != 3
        return false
    end

    ind = model[m4]
    if ind != 4
        return false
    end

    return true
end

"""
Test reaction overloading functions.
"""
function rxn_construction_test(model)
    # note, this test is for iJO1366
    rxn_original = findfirst(model.reactions, "NADH16pp")
    nadh = findfirst(model.metabolites, "nadh_c")
    h_c = findfirst(model.metabolites, "h_c")
    q8 = findfirst(model.metabolites, "q8_c")
    q8h2 = findfirst(model.metabolites, "q8h2_c")
    nad = findfirst(model.metabolites, "nad_c")
    h_p = findfirst(model.metabolites, "h_p")
    
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

"""
Test add!, rm!, and fix_model! functions.
"""
function test_model_manipulations()  
    m1 = Metabolite("m1")
    m2 = Metabolite("m2")
    m3 = Metabolite("m3")
    m4 = Metabolite("m4")
    mets = [m1, m2, m3, m4]
    m5 = Metabolite("m5")
    m6 = Metabolite("m6")
    m7 = Metabolite("m7")
    
    g1 = Gene("g1")
    g2 = Gene("g2")
    g3 = Gene("g3")
    g4 = Gene("g4")
    genes = [g1, g2, g3, g4]
    g5 = Gene("g5")
    g6 = Gene("g6")
    g7 = Gene("g7")
    
    r1 = Reaction("r1", Dict(m1=>-1.0, m2=>1.0) ,"for")
    r2 = Reaction("r2", Dict(m2=>-2.0, m3=>1.0) ,"bidir")
    r2.grr = [[g2], [g1, g3]]
    r3 = Reaction("r3", Dict(m1=>-1.0, m4=>2.0) ,"rev")
    r4 = Reaction("r4", Dict(m1=>-5.0, m4=>2.0) ,"rev")
    r5 = Reaction("r5", Dict(m1=>-11.0, m4=>2.0, m3=>2.0) ,"rev")
    
    rxns = [r1, r2]

    model = Model()
    model.id = "model"
    model.reactions = rxns
    model.metabolites = mets
    model.genes = genes
    
    ### reactions
    add!(model, [r3, r4])
    if length(model.reactions) != 4
        return false
    end
    add!(model, r5)
    if length(model.reactions) != 5
        return false
    end
    rm!(model, [r5, r4])
    if length(model.reactions) != 3
        return false
    end
    rm!(model, r1)
    if length(model.reactions) != 2
        return false
    end

    ### metabolites
    add!(model, [m5, m6])
    if length(model.metabolites) != 6
        return false
    end

    add!(model, m7)
    if length(model.metabolites) != 7
        return false
    end
    rm!(model, [m5, m4])
    if length(model.metabolites) != 5
        return false
    end
    rm!(model, m1)
    if length(model.metabolites) != 4
        return false
    end

    ### genes
    add!(model, [g5, g6])
    if length(model.genes) != 6
        return false
    end

    add!(model, g7)
    if length(model.genes) != 7
        return false
    end
    rm!(model, [g5, g4])
    if length(model.genes) != 5
        return false
    end
    rm!(model, g1)
    if length(model.genes) != 4
        return false
    end
    
    fix_model!(model)
    if length(model.reactions) != 2 && length(model.metabolites) != 4 && length(model.genes) != 3
        return false
    end
    
    return true
end

# function fba_test(model)    
#     biomass_rxn = findfirst(model.rxns, "BIOMASS_Ec_iJO1366_WT_53p95M")
#     solobj = CobraTools.fba(model, biomass_rxn)
#     return solobj.objective ≈ 0.9865144469529787
# end

# function pfba_test(model)    
#     biomass_rxn = findfirst(model.rxns, "BIOMASS_Ec_iJO1366_WT_53p95M")
#     solobj = CobraTools.pfba(model, biomass_rxn)
#     return solobj.objective ≈ 15546.145490407944
# end

# function atom_test(model)
#     biomass_rxn = findfirst(model.rxns, "BIOMASS_Ec_iJO1366_WT_53p95M")
#     pfbasol = CobraTools.pfba(model, biomass_rxn)
#     ad = CobraTools.atom_exchange(pfbasol)
#     return ad["C"]/ad["H"] ≈ 0.6362486422376349
# end