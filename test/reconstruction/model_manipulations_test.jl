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
    model.reactions = OrderedDict(r.id=>r for r in rxns)
    model.metabolites = OrderedDict(m.id=>m for m in mets)
    model.genes = OrderedDict(g.id => g for g in genes)

    ### reactions
    add!(model, [r3, r4])
    @test length(model.reactions) == 4

    add!(model, r5)
    @test length(model.reactions) == 5

    rm!(model, [r5, r4])
    @test length(model.reactions) == 3

    rm!(model, r1)
    @test length(model.reactions) == 2

    ### metabolites
    add!(model, [m5, m6])
    @test length(model.metabolites) == 6

    add!(model, m7)
    @test length(model.metabolites) == 7

    rm!(model, [m5, m4])
    @test length(model.metabolites) == 5

    rm!(model, m1)
    @test length(model.metabolites) == 4

    ### genes
    add!(model, [g5, g6])
    @test length(model.genes) == 6

    add!(model, g7)
    @test length(model.genes) == 7

    rm!(model, [g5, g4])
    @test length(model.genes) == 5

    rm!(model, g1)
    @test length(model.genes) == 4

end
