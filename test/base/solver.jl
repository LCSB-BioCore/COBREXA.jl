
@testset "Solve LP" begin
cp = test_simpleLP()
optimizer = GLPK.Optimizer
(lp, x) = solveLP(cp, optimizer)
@test termination_status(lp) === MOI.OPTIMAL
sol = JuMP.value.(x)
@test sol ≈ [1., 2.]

optimizer = Clp.Optimizer
(lp, x) = solveLP(cp, optimizer)
@test termination_status(lp) === MOI.OPTIMAL
sol = JuMP.value.(x)
@test sol ≈ [1., 2.]
end

