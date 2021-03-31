@testset "LinearModel simple functions" begin
    cp = test_LP()
    @test n_reactions(cp) == 3
    @test n_metabolites(cp) == 4
    @test n_coupling_constraints(cp) == 0

    cp2 = test_LP()
    @test isequal(cp, cp2)
    cp2.S[1] = 1
    @test !isequal(cp, cp2)
    @test isequal(cp, copy(cp))

    cp = test_coupledLP()
    @test n_coupling_constraints(cp) == 2000
    @test isequal(cp, copy(cp))
end
