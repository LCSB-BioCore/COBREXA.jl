
@testset "Solve LP" begin
    cp = test_simpleLP()
    optimizer = Tulip.Optimizer
    lp = optimize_model(cp, optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(lp[:x])
    @test sol ≈ [1.0, 2.0]

    optimizer = Clp.Optimizer
    lp = optimize_model(cp, optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(lp[:x])
    @test sol ≈ [1.0, 2.0]
end
