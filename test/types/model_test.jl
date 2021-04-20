@testset "Model" begin
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
    r1.notes = Dict("notes" => ["blah", "blah"])
    r1.annotation = Dict("sboterm" => "sbo", "biocyc" => ["ads", "asds"])
    r1.objective_coefficient = 1.0

    r2 = Reaction("r2", Dict(m1 => -2.0, m4 => 1.0), :reverse)
    r3 = Reaction("r3", Dict(m3 => -1.0, m4 => 1.0), :forward)
    r4 = Reaction("r4", Dict(m3 => -1.0, m4 => 1.0), :bidirectional)
    r4.annotation = Dict("sboterm" => "sbo", "biocyc" => ["ads", "asds"])

    mets = [m1, m2, m3, m4]
    genes = [g1, g2, g3]
    rxns = [r1, r2, r3, r4]

    model = StandardModel()
    model.id = "model"
    model.reactions = rxns
    model.metabolites = mets
    model.genes = genes

    @test_broken sprint(show, MIME("text/plain"), model) ==
                 "Constraint based model: model\nNumber of reactions: 4\nNumber of metabolites: 4\nNumber of genes: 3\n" # will fix after #92 and PR #93

    @test model[r2] == 2

    @test model[g3] == 3

    @test model[m4] == 4
end
