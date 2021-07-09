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
    genes = [g1, g2, g3, g4]
    g5 = Gene("g5")
    g6 = Gene("g6")
    g7 = Gene("g7")

    r1 = Reaction("r1", Dict(m1.id => -1.0, m2.id => 1.0), :forward)
    r2 = Reaction("r2", Dict(m2.id => -2.0, m3.id => 1.0), :bidirectional)
    r2.grr = [["g2"], ["g1", "g3"]]
    r3 = Reaction("r3", Dict(m1.id => -1.0, m4.id => 2.0), :reverse)
    r4 = Reaction("r4", Dict(m1.id => -5.0, m4.id => 2.0), :reverse)
    r5 = Reaction("r5", Dict(m1.id => -11.0, m4.id => 2.0, m3.id => 2.0), :reverse)

    rxns = [r1, r2]

    model = StandardModel()
    model.id = "model"
    model.reactions = OrderedDict(r.id => r for r in rxns)
    model.metabolites = OrderedDict(m.id => m for m in mets)
    model.genes = OrderedDict(g.id => g for g in genes)

    # change bound tests - in place
    change_bound!(model, "r2"; lower_bound = -10, upper_bound = 10)
    @test model.reactions["r2"].lb == -10
    @test model.reactions["r2"].ub == 10
    change_bound!(model, "r2"; lower_bound = -100)
    @test model.reactions["r2"].lb == -100
    @test model.reactions["r2"].ub == 10
    change_bound!(model, "r2"; upper_bound = 111)
    @test model.reactions["r2"].lb == -100
    @test model.reactions["r2"].ub == 111

    change_bounds!(
        model,
        ["r1", "r2"];
        lower_bounds = [-110, -220],
        upper_bounds = [110.0, 220.0],
    )
    @test model.reactions["r1"].lb == -110
    @test model.reactions["r1"].ub == 110
    @test model.reactions["r2"].lb == -220
    @test model.reactions["r2"].ub == 220

    # change bound - new model
    new_model = change_bound(model, "r2"; lower_bound = -10, upper_bound = 10)
    @test new_model.reactions["r2"].lb == -10
    @test new_model.reactions["r2"].ub == 10

    new_model = change_bound(model, "r2"; lower_bound = -10)
    @test new_model.reactions["r2"].lb == -10
    @test new_model.reactions["r2"].ub == 220

    new_model = change_bounds(
        model,
        ["r1", "r2"];
        lower_bounds = [-10, -20],
        upper_bounds = [10.0, 20.0],
    )
    @test new_model.reactions["r1"].lb == -10
    @test new_model.reactions["r1"].ub == 10
    @test new_model.reactions["r2"].lb == -20
    @test new_model.reactions["r2"].ub == 20

    new_model = change_bounds(
        model,
        ["r1", "r2"];
        lower_bounds = [-10, nothing],
        upper_bounds = [nothing, 20.0],
    )
    @test new_model.reactions["r1"].lb == -10
    @test new_model.reactions["r1"].ub == 110
    @test new_model.reactions["r2"].lb == -220
    @test new_model.reactions["r2"].ub == 20

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
end
