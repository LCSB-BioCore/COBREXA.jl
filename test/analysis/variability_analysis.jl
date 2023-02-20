@testset "Flux variability analysis" begin
    cp = test_simpleLP()
    optimizer = Tulip.Optimizer
    fluxes = flux_variability_analysis(cp, optimizer) |> result

    @test size(fluxes) == (2, 2)
    @test fluxes ≈ [
        1.0 1.0
        2.0 2.0
    ]

    rates = variability_analysis(cp, optimizer) |> result
    @test fluxes == rates

    fluxes = flux_variability_analysis(cp, optimizer, reaction_indexes = [2]) |> result

    @test size(fluxes) == (1, 2)
    @test isapprox(fluxes, [2 2], atol = TEST_TOLERANCE)

    # a special testcase for slightly sub-optimal FVA (gamma<1)
    cp = MatrixModel(
        [-1.0 -1.0 -1.0],
        [0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, -1.0],
        1.0 * ones(3),
        ["r$x" for x = 1:3],
        ["m1"],
    )
    fluxes = flux_variability_analysis(cp, optimizer) |> result
    @test isapprox(
        fluxes,
        [
            1.0 1.0
            0.0 0.0
            -1.0 -1.0
        ],
        atol = TEST_TOLERANCE,
    )
    fluxes = flux_variability_analysis(cp, optimizer; bounds = gamma_bounds(0.5)) |> result
    @test isapprox(
        fluxes,
        [
            0.5 1.0
            0.0 0.5
            -1.0 -0.5
        ],
        atol = TEST_TOLERANCE,
    )
    fluxes = flux_variability_analysis(cp, optimizer; bounds = _ -> (0, Inf)) |> result
    @test isapprox(
        fluxes,
        [
            0.0 1.0
            0.0 1.0
            -1.0 0.0
        ],
        atol = TEST_TOLERANCE,
    )

    @test isempty(
        flux_variability_analysis(cp, Tulip.Optimizer, reaction_ids = String[]) |> result,
    )
    @test_throws DomainError flux_variability_analysis(
        cp,
        Tulip.Optimizer,
        reaction_indexes = [-1],
    )
    @test_throws DomainError flux_variability_analysis(
        cp,
        Tulip.Optimizer,
        reaction_ids = ["not a reaction!"],
    )
end

@testset "Parallel FVA" begin
    cp = test_simpleLP()

    fluxes =
        flux_variability_analysis(
            cp,
            Tulip.Optimizer;
            workers = W,
            reaction_indexes = [1, 2],
        ) |> result
    @test isapprox(
        fluxes,
        [
            1.0 1.0
            2.0 2.0
        ],
        atol = TEST_TOLERANCE,
    )
end

@testset "Flux variability analysis with ObjectModel" begin
    model = load_model(ObjectModel, model_paths["e_coli_core.json"])
    mins, maxs =
        flux_variability_analysis_dict(
            model,
            Tulip.Optimizer;
            bounds = objective_bounds(0.99),
            modifications = [
                modify_optimizer_attribute("IPM_IterationsLimit", 500),
                modify_constraint("EX_glc__D_e"; lower_bound = -10, upper_bound = -10),
                modify_constraint("EX_o2_e"; lower_bound = 0.0, upper_bound = 0.0),
            ],
        ) |> result

    @test isapprox(maxs["EX_ac_e"]["EX_ac_e"], 8.5185494, atol = TEST_TOLERANCE)
    @test isapprox(mins["EX_ac_e"]["EX_ac_e"], 7.4483887, atol = TEST_TOLERANCE)
end
