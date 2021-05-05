@testset "Flux variability analysis" begin
    cp = test_simpleLP()
    optimizer = GLPK.Optimizer
    fluxes = flux_variability_analysis(cp, optimizer)

    @test size(fluxes) == (2, 2)
    @test fluxes ≈ [
        1.0 1.0
        2.0 2.0
    ]

    fluxes = flux_variability_analysis(cp, [2], optimizer)

    @test size(fluxes) == (1, 2)
    @test fluxes == Matrix{Float64}([2 2])

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
    @test fluxes ≈ [
        1.0 1.0
        0.0 0.0
        -1.0 -1.0
    ]
    fluxes = flux_variability_analysis(cp, optimizer; bounds = gamma_bounds(0.5))
    @test fluxes ≈ [
        0.5 1.0
        0.0 0.5
        -1.0 -0.5
    ]
    fluxes = flux_variability_analysis(cp, optimizer; bounds = _ -> (0, Inf))
    @test fluxes ≈ [
        0.0 1.0
        0.0 1.0
        -1.0 0.0
    ]

    @test isempty(flux_variability_analysis(cp, Vector{Int}(), GLPK.Optimizer))
    @test_throws DomainError flux_variability_analysis(cp, [-1], GLPK.Optimizer)
    @test_throws DomainError flux_variability_analysis(cp, [99999999], GLPK.Optimizer)
end

@testset "Parallel FVA" begin
    cp = test_simpleLP()
    pids = addprocs(2, topology = :master_worker)
    @everywhere using COBREXA, GLPK
    fluxes = flux_variability_analysis(cp, [1, 2], GLPK.Optimizer; workers = pids)
    @test fluxes ≈ [
        1.0 1.0
        2.0 2.0
    ]
    rmprocs(pids)
end

@testset "Flux variability analysis with StandardModel" begin
    model_path = download_data_file(
        "http://bigg.ucsd.edu/static/models/e_coli_core.json",
        joinpath("data", "e_coli_core.json"),
        "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
    )

    model = load_model(StandardModel, model_path)
    mins, maxs = flux_variability_analysis_dict(
        model,
        Tulip.Optimizer;
        bounds = objective_bounds(0.99),
        modifications = [
            change_solver_attribute("IPM_IterationsLimit", 500),
            change_constraint("EX_glc__D_e", -10, -10),
            change_constraint("EX_o2_e", 0.0, 0.0),
        ],
    )

    @test isapprox(maxs["EX_ac_e"]["EX_ac_e"], 8.518549434876208, atol = TEST_TOLERANCE)
    @test isapprox(mins["EX_ac_e"]["EX_ac_e"], 7.448388738973361, atol = TEST_TOLERANCE)
end
