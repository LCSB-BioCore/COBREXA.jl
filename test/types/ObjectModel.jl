@testset "ObjectModel generic interface" begin
    # create a small model
    m1 = Metabolite("m1")
    m1.formula = "C2H3"
    m1.compartment = "cytosol"
    m2 = Metabolite("m2")
    m2.formula = "H3C2"
    m3 = Metabolite("m3")
    m3.charge = -1
    m4 = Metabolite("m4")
    m4.notes = Dict("confidence" => ["iffy"])
    m4.annotations = Dict("sbo" => ["blah"])

    g1 = Gene("g1")
    g2 = Gene("g2")
    g2.notes = Dict("confidence" => ["iffy"])
    g2.annotations = Dict("sbo" => ["blah"])
    g3 = Gene("g3")

    r1 = Reaction(id = "r1")
    r1.metabolites = Dict(m1.id => -1.0, m2.id => 1.0)
    r1.lower_bound = -100.0
    r1.upper_bound = 100.0
    r1.gene_associations = [
        Isozyme(stoichiometry = Dict("g1" => 1, "g2" => 1)),
        Isozyme(stoichiometry = Dict("g3" => 1)),
    ]
    r1.subsystem = "glycolysis"
    r1.notes = Dict("notes" => ["blah", "blah"])
    r1.annotations = Dict("sboterm" => ["sbo"], "biocyc" => ["ads", "asds"])

    r2 = Reaction("r2", Dict(m1.id => -2.0, m4.id => 1.0), :reverse)
    r3 = Reaction("r3", Dict(m3.id => -1.0, m4.id => 1.0), :forward)
    r4 = Reaction("r4", Dict(m3.id => -1.0, m4.id => 1.0), :bidirectional)
    r4.annotations = Dict("sboterm" => ["sbo"], "biocyc" => ["ads", "asds"])

    mets = [m1, m2, m3, m4]
    gs = [g1, g2, g3]
    rxns = [r1, r2, r3, r4]

    model = ObjectModel(id = "model")
    model.reactions = OrderedDict(r.id => r for r in rxns)
    model.metabolites = OrderedDict(m.id => m for m in mets)
    model.genes = OrderedDict(g.id => g for g in gs)
    model.objective = Dict("r1" => 1.0)

    @test contains(sprint(show, MIME("text/plain"), model), "ObjectModel")

    @test "r1" in reactions(model)
    @test "m4" in metabolites(model)
    @test "g2" in genes(model)
    @test n_variables(model) == 4
    @test n_metabolites(model) == 4
    @test n_genes(model) == 3

    S_test = spzeros(4, 4)
    S_test[1, 1] = -1.0
    S_test[2, 1] = 1.0
    S_test[1, 2] = -2.0
    S_test[4, 2] = 1.0
    S_test[3, 3] = -1.0
    S_test[4, 3] = 1.0
    S_test[3, 4] = -1.0
    S_test[4, 4] = 1.0
    @test S_test == stoichiometry(model)

    lb_test = spzeros(4)
    lb_test[1] = -100.0
    lb_test[2] = -1000.0
    lb_test[3] = 0.0
    lb_test[4] = -1000.0
    ub_test = spzeros(4)
    ub_test[1] = 100.0
    ub_test[2] = 0.0
    ub_test[3] = 1000.0
    ub_test[4] = 1000.0
    lbs, ubs = bounds(model)
    @test lb_test == lbs
    @test ub_test == ubs

    @test balance(model) == spzeros(n_metabolites(model))

    obj_test = spzeros(4)
    obj_test[1] = 1.0
    @test objective(model) == obj_test

    @test all(
        occursin.(
            ["g1", "g2", "g3"],
            Ref(
                COBREXA.Internal.unparse_grr(
                    String,
                    reaction_gene_association(model, "r1"),
                ),
            ),
        ),
    )
    @test isnothing(reaction_gene_association(model, "r2"))

    @test metabolite_formula(model, "m2")["C"] == 2
    @test isnothing(metabolite_formula(model, "m3"))

    @test metabolite_charge(model, "m3") == -1
    @test isnothing(metabolite_charge(model, "m2"))

    @test metabolite_compartment(model, "m1") == "cytosol"
    @test isnothing(metabolite_compartment(model, "m2"))

    @test reaction_subsystem(model, "r1") == "glycolysis"
    @test isnothing(reaction_subsystem(model, "r2"))

    @test metabolite_notes(model, "m4")["confidence"] == ["iffy"]
    @test metabolite_annotations(model, "m4")["sbo"] == ["blah"]
    @test isempty(metabolite_notes(model, "m3"))
    @test isempty(metabolite_annotations(model, "m3"))

    @test gene_notes(model, "g2")["confidence"] == ["iffy"]
    @test gene_annotations(model, "g2")["sbo"] == ["blah"]
    @test isempty(gene_notes(model, "g1"))
    @test isempty(gene_annotations(model, "g1"))

    @test reaction_notes(model, "r1")["notes"] == ["blah", "blah"]
    @test reaction_annotations(model, "r1")["biocyc"] == ["ads", "asds"]
    @test isempty(reaction_notes(model, "r2"))
    @test isempty(reaction_annotations(model, "r2"))

    @test reaction_stoichiometry(model, "r1") == Dict("m1" => -1.0, "m2" => 1.0)

    # To do: test convert
    same_model = convert(ObjectModel, model)
    @test same_model == model

    jsonmodel = convert(JSONModel, model)
    stdmodel = convert(ObjectModel, jsonmodel)
    @test issetequal(reactions(jsonmodel), reactions(stdmodel))
    @test issetequal(genes(jsonmodel), genes(stdmodel))
    @test issetequal(metabolites(jsonmodel), metabolites(stdmodel))
    jlbs, jubs = bounds(jsonmodel)
    slbs, subs = bounds(stdmodel)
    @test issetequal(jlbs, slbs)
    @test issetequal(jubs, subs)
    jS = stoichiometry(jsonmodel)
    sS = stoichiometry(stdmodel)
    j_r1_index = findfirst(x -> x == "r1", reactions(jsonmodel))
    s_r1_index = findfirst(x -> x == "r1", reactions(stdmodel))
    j_m1_index = findfirst(x -> x == "m1", metabolites(jsonmodel))
    j_m2_index = findfirst(x -> x == "m2", metabolites(jsonmodel))
    s_m1_index = findfirst(x -> x == "m1", metabolites(stdmodel))
    s_m2_index = findfirst(x -> x == "m2", metabolites(stdmodel))
    @test jS[j_m1_index, j_r1_index] == sS[s_m1_index, s_r1_index]
    @test jS[j_m2_index, j_r1_index] == sS[s_m2_index, s_r1_index]
end
