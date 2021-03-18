@testset "Flux variability analysis" begin
    cp = test_simpleLP()
    optimizer = GLPK.Optimizer
    fluxes = fluxVariabilityAnalysis(cp, optimizer)

    @test size(fluxes) == (2, 2)
    @test fluxes ≈ [
        1.0 1.0
        2.0 2.0
    ]

    fluxes = fluxVariabilityAnalysis(cp, [2], optimizer)

    @test size(fluxes) == (1, 2)
    @test fluxes == Array{Float64,2}([2 2])

    cp = LinearModel(
        [
            -1.0 -1.0 -1.0
        ],
        [0.],
        [1.0, 0., 0.],
        [0., 0., -1.0],
        1.0 * ones(3),
        ["r$x" for x = 1:3],
        ["m1"],
    )
    fluxes = fluxVariabilityAnalysis(cp, optimizer)
    @test fluxes ≈ [
        1.0   1.0;
        0.0   0.0;
        -1.0  -1.0;
    ]
    fluxes = fluxVariabilityAnalysis(cp, optimizer, 0.5)
    @test fluxes ≈ [
        0.5   1.0;
        0.0   0.5;
        -1.0  -0.5;
    ]
    fluxes = fluxVariabilityAnalysis(cp, optimizer, 0.)
    @test fluxes ≈ [
        0.   1.0;
        0.0   1.0;
        -1.0  0.;
    ]
end

@testset "Parallel FVA" begin
    cp = test_simpleLP()
    pids = addprocs(2, topology = :master_worker)
    @everywhere using COBREXA, GLPK
    fluxes = parFVA(cp, [1, 2], GLPK.Optimizer, pids)
    @test fluxes ≈ [
        1.0 1.0
        2.0 2.0
    ]
    @test_throws ArgumentError parFVA(cp, [99999999], GLPK.Optimizer, pids)
    rmprocs(pids)
end
