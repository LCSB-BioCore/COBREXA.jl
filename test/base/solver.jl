
@testset "Solve LP" begin
    cp = test_simpleLP()
    optimizer = GLPK.Optimizer
    (lp, x) = solveLP(cp, optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(x)
    @test sol ≈ [1.0, 2.0]

    optimizer = Clp.Optimizer
    (lp, x) = solveLP(cp, optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(x)
    @test sol ≈ [1.0, 2.0]
end
