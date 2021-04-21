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
    r1.metabolites = Dict(m1 => -1.0, m2 => 1.0)
    r1.lb = -100.0
    r1.ub = 100.0
    r1.grr = [[g1, g2], [g3]]
    r1.subsystem = "glycolysis"
    r1.notes = Dict("notes" => ["blah", "blah"])
    r1.annotation = Dict("sboterm" => "sbo", "biocyc" => ["ads", "asds"])
    r1.objective_coefficient = 1.0

    @test sprint(show, MIME("text/plain"), r1) ==
          "Reaction ID: r1\nName: reaction 1\nReaction equation: 1.0 m1 ⟷  1.0 m2\nLower bound: -100.0\nUpper bound: 100.0\nSubsystem: glycolysis\nGene reaction rule: (g1 and g2) or (g3)\nNotes: \n\tnotes: blah, blah\nAnnotation: \n\tsboterm: sbo\n\tbiocyc: ads, asds\nFields: id, name, metabolites, lb, ub, grr, subsystem, notes, annotation, objective_coefficient\n"

    rlongfor = Reaction(
        "rlongfor",
        Dict(
            m1 => -1.0,
            m2 => -1.0,
            m3 => -1.0,
            m4 => -1.0,
            m5 => -1.0,
            m6 => -1.0,
            m7 => 1.0,
            m8 => 1.0,
            m9 => 1.0,
            m10 => 1.0,
            m11 => 1.0,
            m12 => 1.0,
        ),
        :forward,
    )
    @test occursin("...", sprint(show, MIME("text/plain"), rlongfor)) # from dictionaries so order is not consistent.

    rlongrev = Reaction(
        "rlongrev",
        Dict(
            m1 => -1.0,
            m2 => -1.0,
            m3 => -1.0,
            m4 => -1.0,
            m5 => -1.0,
            m6 => -1.0,
            m7 => 1.0,
            m8 => 1.0,
            m9 => 1.0,
            m10 => 1.0,
            m11 => 1.0,
            m12 => 1.0,
        ),
        :reverse,
    )
    @test occursin("...", sprint(show, MIME("text/plain"), rlongrev))

    r2 = Reaction("r2", Dict(m1 => -2.0, m4 => 1.0), :reverse)
    @test r2.lb == -1000.0 && r2.ub == 0.0

    r3 = Reaction("r3", Dict(m3 => -1.0, m4 => 1.0), :forward)
    @test r3.lb == 0.0 && r3.ub == 1000.0

    rxns = [r1, r2, r3]

    @test sprint(show, MIME("text/plain"), rxns) ==
          "Reaction vector of length: : 3\nEach reaction has fields: id, name, metabolites, lb, ub, grr, subsystem, notes, annotation, objective_coefficient\n"

    @test rxns[r3] == 3

    rr = findfirst(rxns, "r2")
    @test rr.id == r2.id

    r4 = Reaction("r4", Dict(m3 => -1.0, m4 => 1.0), :bidirectional)
    r4.annotation = Dict("sboterm" => "sbo", "biocyc" => ["ads", "asds"])
    @test r4.lb == -1000.0 && r4.ub == 1000.0

    dup, ind = check_duplicate_annotations(rxns, r4)
    @test dup && ind == 1

    dup, ind = check_duplicate_reaction(rxns, r4)
    @test dup && ind == 3

    dup, ind = check_duplicate_annotations(rxns, r2)
    @test !dup && ind == -1

    r5 = Reaction("r5", Dict(m3 => -11.0, m4 => 1.0), :bidirectional)
    dup, ind = check_duplicate_reaction(rxns, r5)
    @test !dup && ind == -1

    bal, d = is_mass_balanced(r1)
    @test bal

    @test isnothing(findfirst(rxns, "nope"))
end
