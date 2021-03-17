@testset "Flux balance analysis" begin
    cp = test_simpleLP()
    optimizer = GLPK.Optimizer
    (lp, x) = fluxBalanceAnalysis(cp, optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(x)
    @test sol ≈ [1., 2.]

    optimizer = Clp.Optimizer
    (lp, x) = fluxBalanceAnalysis(cp, optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(x)
    @test sol ≈ [1., 2.]

    # test the maximization of the objective
    cp = test_simpleLP2()
    optimizer = GLPK.Optimizer
    (lp, x) = fluxBalanceAnalysis(cp, optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(x)
    @test sol ≈ [-1., 2.]

    # test with a more biologically meaningfull model
    modelPath = joinpath("data", "fba.mat")
    cp = loadModel(modelPath, "model")
    file = matopen(modelPath)
    expectedOptimum = matread(modelPath)["optimum"]
    close(file)

    optimizer = GLPK.Optimizer
    (lp, x) = fluxBalanceAnalysis(cp, optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(x)
    @test objective_value(lp) ≈ expectedOptimum
    @test cp.c' * sol ≈ expectedOptimum
end