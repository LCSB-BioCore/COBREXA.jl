@testset "Flux balance analysis" begin
    cp = test_simpleLP()
    (lp, x) = fluxBalanceAnalysis(cp, GLPK.Optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(x)
    @test sol ≈ [1.0, 2.0]

    (lp, x) = fluxBalanceAnalysis(cp, Clp.Optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(x)
    @test sol ≈ [1.0, 2.0]

    # test the maximization of the objective
    cp = test_simpleLP2()
    (lp, x) = fluxBalanceAnalysis(cp, GLPK.Optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(x)
    @test sol ≈ [-1.0, 2.0]

    # test with a more biologically meaningfull model
    modelPath = joinpath("data", "fba.mat")
    cp = loadModel(modelPath, "model")
    file = matopen(modelPath)
    expectedOptimum = matread(modelPath)["optimum"]
    close(file)

    (lp, x) = fluxBalanceAnalysis(cp, GLPK.Optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(x)
    @test objective_value(lp) ≈ expectedOptimum
    @test cp.c' * sol ≈ expectedOptimum

    # test the "nicer output" variants
    @test_broken false # reminder to implement these methods
    # fluxesVec = fluxBalanceAnalysisVec(cp, GLPK.Optimizer)
    # @test_broken all(fluxesVec .== sol)
    # fluxesDict = fluxBalanceAnalysisDict(cp, GLPK.Optimizer)
    # rxns = reactions(cp)
    # @test all([fluxesDict[rxns[i]] == sol[i] for i in eachindex(rxns)])
end
