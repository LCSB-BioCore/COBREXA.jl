
@testset "Solve LP" begin
    cp = test_simpleLP()
    lp = optimize_model(cp, Tulip.Optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(lp[:x])
    @test sol ≈ [1.0, 2.0]
end
