@testset "Reaction" begin
    m1 = Metabolite("m1")
    m1.formula = "C2H3"
    m2 = Metabolite("m2")
    m2.formula = "H3C2"
    m3 = Metabolite("m3")
    m4 = Metabolite("m4")
    m5 = Metabolite("m5")
    m6 = Metabolite("m6")
    m7 = Metabolite("m7")
    m8 = Metabolite("m8")
    m9 = Metabolite("m9")
    m10 = Metabolite("m10")
    m11 = Metabolite("m11")
    m12 = Metabolite("m12")

    g1 = Gene("g1")
    g2 = Gene("g2")
    g3 = Gene("g3")

    r1 = Reaction()
    r1.id = "r1"
    r1.name = "reaction 1"
    r1.metabolites = Dict(m1.id => -1.0, m2.id => 1.0)
    r1.lb = -100.0
    r1.ub = 100.0
    r1.grr = [["g1", "g2"], ["g3"]]
    r1.subsystem = "glycolysis"
    r1.notes = Dict("notes" => ["blah", "blah"])
    r1.annotation = Dict("sboterm" => ["sbo"], "biocyc" => ["ads", "asds"])
    r1.objective_coefficient = 1.0

    @test sprint(show, MIME("text/plain"), r1) ==
          "Reaction ID: r1\nReaction name: reaction 1\nReaction subsystem: glycolysis\n1.0 m1 ⟷  1.0 m2\nLower bound: -100.0\nUpper bound: 100.0\nGenes: (g1 and g2) or (g3)\nE.C. number: \n"

    rlongfor = Reaction(
        "rlongfor",
        Dict(
            m1.id => -1.0,
            m2.id => -1.0,
            m3.id => -1.0,
            m4.id => -1.0,
            m5.id => -1.0,
            m6.id => -1.0,
            m7.id => 1.0,
            m8.id => 1.0,
            m9.id => 1.0,
            m10.id => 1.0,
            m11.id => 1.0,
            m12.id => 1.0,
        ),
        :forward,
    )
    @test occursin("...", sprint(show, MIME("text/plain"), rlongfor)) # from dictionaries so order is not consistent.

    rlongrev = Reaction(
        "rlongrev",
        Dict(
            m1.id => -1.0,
            m2.id => -1.0,
            m3.id => -1.0,
            m4.id => -1.0,
            m5.id => -1.0,
            m6.id => -1.0,
            m7.id => 1.0,
            m8.id => 1.0,
            m9.id => 1.0,
            m10.id => 1.0,
            m11.id => 1.0,
            m12.id => 1.0,
        ),
        :reverse,
    )
    @test occursin("...", sprint(show, MIME("text/plain"), rlongrev))

    r2 = Reaction("r2", Dict(m1.id => -2.0, m4.id => 1.0), :reverse)
    @test r2.lb == -1000.0 && r2.ub == 0.0

    r3 = Reaction("r3", Dict(m3.id => -1.0, m4.id => 1.0), :forward)
    @test r3.lb == 0.0 && r3.ub == 1000.0

    rxns = [r1, r2, r3]
    rd = OrderedDict(r.id => r for r in rxns)
    @test sprint(show, MIME("text/plain"), rxns) == "Reaction set of length: 3\n"

    r4 = Reaction("r4", Dict(m3.id => -1.0, m4.id => 1.0), :bidirectional)
    r4.annotation = Dict("sboterm" => ["sbo"], "biocyc" => ["ads", "asds"])
    @test r4.lb == -1000.0 && r4.ub == 1000.0

    dup, ind = check_duplicate_annotations(r4, rd)
    @test dup && ind == "r1"

    dup, ind = check_duplicate_reaction(r4, rd)
    @test dup && ind == "r3"

    dup, ind = check_duplicate_annotations(r2, rd)
    @test !dup && ind == ""

    r5 = Reaction("r5", Dict(m3.id => -11.0, m4.id => 1.0), :bidirectional)
    dup, ind = check_duplicate_reaction(r5, rd)
    @test !dup && ind == ""
end
