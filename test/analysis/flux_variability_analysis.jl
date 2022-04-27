@testset "Flux variability analysis" begin
    cp = test_simpleLP()
    optimizer = Tulip.Optimizer
    fluxes = flux_variability_analysis(cp, optimizer)

    @test size(fluxes) == (2, 2)
    @test fluxes ≈ [
        1.0 1.0
        2.0 2.0
    ]

    rates = reaction_variability_analysis(cp, optimizer)
    @test fluxes == rates

    fluxes = flux_variability_analysis(cp, [2], optimizer)

    @test size(fluxes) == (1, 2)
    @test isapprox(fluxes, [2 2], atol = TEST_TOLERANCE)

    # a special testcase for slightly sub-optimal FVA (gamma<1)
    cp = CoreModel(
        [-1.0 -1.0 -1.0],
        [0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, -1.0],
        1.0 * ones(3),
        ["r$x" for x = 1:3],
        ["m1"],
    )
    fluxes = flux_variability_analysis(cp, optimizer)
    @test isapprox(
        fluxes,
        [
            1.0 1.0
            0.0 0.0
            -1.0 -1.0
        ],
        atol = TEST_TOLERANCE,
    )
    fluxes = flux_variability_analysis(cp, optimizer; bounds = gamma_bounds(0.5))
    @test isapprox(
        fluxes,
        [
            0.5 1.0
            0.0 0.5
            -1.0 -0.5
        ],
        atol = TEST_TOLERANCE,
    )
    fluxes = flux_variability_analysis(cp, optimizer; bounds = _ -> (0, Inf))
    @test isapprox(
        fluxes,
        [
            0.0 1.0
            0.0 1.0
            -1.0 0.0
        ],
        atol = TEST_TOLERANCE,
    )

    @test isempty(flux_variability_analysis(cp, Vector{Int}(), Tulip.Optimizer))
    @test_throws DomainError flux_variability_analysis(cp, [-1], Tulip.Optimizer)
    @test_throws DomainError flux_variability_analysis(cp, [99999999], Tulip.Optimizer)
end

@testset "Parallel FVA" begin
    cp = test_simpleLP()

    fluxes = flux_variability_analysis(cp, [1, 2], Tulip.Optimizer; workers = W)
    @test isapprox(
        fluxes,
        [
            1.0 1.0
            2.0 2.0
        ],
        atol = TEST_TOLERANCE,
    )
end

@testset "Flux variability analysis with StandardModel" begin
    model = load_model(StandardModel, model_paths["e_coli_core.json"])
    mins, maxs = flux_variability_analysis_dict(
        model,
        Tulip.Optimizer;
        bounds = objective_bounds(0.99),
        modifications = [
            change_optimizer_attribute("IPM_IterationsLimit", 500),
            change_constraint("EX_glc__D_e"; lb = -10, ub = -10),
            change_constraint("EX_o2_e"; lb = 0.0, ub = 0.0),
        ],
    )

    @test isapprox(maxs["EX_ac_e"]["EX_ac_e"], 8.5185494, atol = TEST_TOLERANCE)
    @test isapprox(mins["EX_ac_e"]["EX_ac_e"], 7.4483887, atol = TEST_TOLERANCE)
end
