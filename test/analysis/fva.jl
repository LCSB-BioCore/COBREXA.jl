@testset "Flux variability analysis" begin
    cp = test_simpleLP()
    optimizer = GLPK.Optimizer
    fluxes = fluxVariabilityAnalysis(cp, optimizer)

    @test size(fluxes) == (2, 2)
    @test fluxes ≈ [1. 1.;
                    2. 2.]

    fluxes = fluxVariabilityAnalysis(cp, [2], optimizer)

    @test size(fluxes) == (1, 2)
    @test fluxes == Array{Float64, 2}([2 2])
end

@testset "Parallel FVA" begin
    cp = test_simpleLP()
    pids = addprocs(2, topology=:master_worker)
    @everywhere using COBREXA, GLPK
    fluxes = parFVA(cp, [1,2], GLPK.Optimizer, pids)
    @test fluxes ≈ [1. 1.;
                    2. 2.]
    @test_throws ArgumentError parFVA(cp, [99999999], GLPK.Optimizer, pids)
    rmprocs(pids)
end
