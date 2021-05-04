@testset "single_knockout" begin
    optimizer = Tulip.Optimizer
    m = StandardModel()
    add!(m, Metabolite("A"))
    add!(m, Metabolite("B"))

    add!(m, Gene("g1"))
    add!(m, Gene("g2"))
    add!(m, Reaction("v1", metabolites = Dict("A" => -1.0, "B" => 1.0), grr = [["g1"]]))
    add!(
        m,
        Reaction("v2", metabolites = Dict("A" => -1.0, "B" => 1.0), grr = [["g1", "g2"]]),
    )
    add!(
        m,
        Reaction("v3", metabolites = Dict("A" => -1.0, "B" => 1.0), grr = [["g1"], ["g2"]]),
    )
    add!(
        m,
        Reaction(
            "v4",
            metabolites = Dict("A" => -1.0, "B" => 1.0),
            grr = [["g1", "g2"], ["g2"]],
        ),
    )

    opt_model = make_optimization_model(m, optimizer)
    knockout("g1")(m, opt_model)

    # Knockout should remove v1
    @test normalized_rhs(opt_model[:lbs][1]) == 0
    @test normalized_rhs(opt_model[:ubs][1]) == 0

    # Knockout should remove [g1, g2] (AND) and thus remove reaction
    @test normalized_rhs(opt_model[:lbs][2]) == 0
    @test normalized_rhs(opt_model[:ubs][2]) == 0

    # Knockout should remove [g1], but keep reaction (OR)
    @test normalized_rhs(opt_model[:lbs][3]) == 1000
    @test normalized_rhs(opt_model[:ubs][3]) == 1000

    # Knockout should remove [g1, g2] (AND), but keep reaction (OR)
    @test normalized_rhs(opt_model[:lbs][4]) == 1000
    @test normalized_rhs(opt_model[:ubs][4]) == 1000
end

@testset "multiple_knockouts" begin
    optimizer = Tulip.Optimizer
    m = StandardModel()
    add!(m, Metabolite("A"))
    add!(m, Metabolite("B"))
    add!(m, Gene("g1"))
    add!(m, Gene("g2"))
    add!(m, Gene("g3"))
    add!(
        m,
        Reaction("v1", metabolites = Dict("A" => -1.0, "B" => 1.0), grr = [["g1"], ["g3"]]),
    )
    add!(
        m,
        Reaction(
            "v2",
            metabolites = Dict("A" => -1.0, "B" => 1.0),
            grr = [["g1", "g2"], ["g3"]],
        ),
    )
    add!(
        m,
        Reaction(
            "v3",
            metabolites = Dict("A" => -1.0, "B" => 1.0),
            grr = [["g1"], ["g2"], ["g3"]],
        ),
    )

    opt_model = make_optimization_model(m, optimizer)
    knockout(["g1", "g3"])(m, opt_model)

    # Reaction 1 should be knocked out, because both
    # gene1 and gene 3 are knocked out
    @test normalized_rhs(opt_model[:lbs][1]) == 0
    @test normalized_rhs(opt_model[:ubs][1]) == 0

    # Reaction 2 should be knocked out, because both
    # [g1, g2] is an AND relationship
    @test normalized_rhs(opt_model[:lbs][1]) == 0
    @test normalized_rhs(opt_model[:ubs][1]) == 0

    # Reaction 3 should stay, because gene2 is still
    # available (the arrays have an OR relationship)
    @test normalized_rhs(opt_model[:lbs][3]) == 1000
    @test normalized_rhs(opt_model[:ubs][3]) == 1000
end
