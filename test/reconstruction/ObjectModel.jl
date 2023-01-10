@testset "Model manipulation" begin
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
    gene_vec = [g1, g2, g3, g4]
    g5 = Gene("g5")
    g6 = Gene("g6")
    g7 = Gene("g7")

    r1 = Reaction("r1", Dict(m1.id => -1.0, m2.id => 1.0), :forward)
    r2 = Reaction("r2", Dict(m2.id => -2.0, m3.id => 1.0), :bidirectional)
    r2.gene_associations = [Isozyme(x) for x in [["g2"], ["g1", "g3"]]]
    r3 = Reaction("r3", Dict(m1.id => -1.0, m4.id => 2.0), :backward)
    r4 = Reaction("r4", Dict(m1.id => -5.0, m4.id => 2.0), :backward)
    r5 = Reaction("r5", Dict(m1.id => -11.0, m4.id => 2.0, m3.id => 2.0), :backward)

    rxns = [r1, r2]

    model = ObjectModel(id = "model")
    model.reactions = OrderedDict(r.id => r for r in rxns)
    model.metabolites = OrderedDict(m.id => m for m in mets)
    model.genes = OrderedDict(g.id => g for g in gene_vec)

    # change bound tests - in place
    change_bound!(model, "r2"; lower = -10, upper = 10)
    @test model.reactions["r2"].lower_bound == -10
    @test model.reactions["r2"].upper_bound == 10
    change_bound!(model, "r2"; lower = -100)
    @test model.reactions["r2"].lower_bound == -100
    @test model.reactions["r2"].upper_bound == 10
    change_bound!(model, "r2"; upper = 111)
    @test model.reactions["r2"].lower_bound == -100
    @test model.reactions["r2"].upper_bound == 111

    change_bounds!(model, ["r1", "r2"]; lower = [-110, -220], upper = [110.0, 220.0])
    @test model.reactions["r1"].lower_bound == -110
    @test model.reactions["r1"].upper_bound == 110
    @test model.reactions["r2"].lower_bound == -220
    @test model.reactions["r2"].upper_bound == 220

    # change bound - new model
    new_model = change_bound(model, "r2"; lower = -10, upper = 10)
    @test new_model.reactions["r2"].lower_bound == -10
    @test new_model.reactions["r2"].upper_bound == 10

    new_model = change_bound(model, "r2"; lower = -10)
    @test new_model.reactions["r2"].lower_bound == -10
    @test new_model.reactions["r2"].upper_bound == 220

    new_model = change_bounds(model, ["r1", "r2"]; lower = [-10, -20], upper = [10.0, 20.0])
    @test new_model.reactions["r1"].lower_bound == -10
    @test new_model.reactions["r1"].upper_bound == 10
    @test new_model.reactions["r2"].lower_bound == -20
    @test new_model.reactions["r2"].upper_bound == 20

    new_model =
        change_bounds(model, ["r1", "r2"]; lower = [-10, nothing], upper = [nothing, 20.0])
    @test new_model.reactions["r1"].lower_bound == -10
    @test new_model.reactions["r1"].upper_bound == 110
    @test new_model.reactions["r2"].lower_bound == -220
    @test new_model.reactions["r2"].upper_bound == 20

    ### reactions
    add_reactions!(model, [r3, r4])
    @test length(model.reactions) == 4

    add_reaction!(model, r5)
    @test length(model.reactions) == 5

    remove_reactions!(model, ["r5", "r4"])
    @test length(model.reactions) == 3

    remove_reaction!(model, "r1")
    @test length(model.reactions) == 2

    ### metabolites
    add_metabolites!(model, [m5, m6])
    @test length(model.metabolites) == 6

    add_metabolite!(model, m7)
    @test length(model.metabolites) == 7

    remove_metabolites!(model, ["m5", "m4"])
    @test length(model.metabolites) == 5

    remove_metabolite!(model, "m1")
    @test length(model.metabolites) == 4

    ### genes
    add_genes!(model, [g5, g6])
    @test length(model.genes) == 6

    add_gene!(model, g7)
    @test length(model.genes) == 7

    remove_genes!(model, ["g5", "g4"])
    @test length(model.genes) == 5

    remove_gene!(model, "g1")
    @test length(model.genes) == 4

    ### objective
    change_objective!(model, "r2")
    @test model.objective["r2"] == 1.0
end
