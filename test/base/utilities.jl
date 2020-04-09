@testset "LinearModel simple functions" begin
    cp = test_LP()
    @test nReactions(cp) == 3
    @test nMetabolites(cp) == 4
    @test nCouplingConstraints(cp) == 0

    cp2 = test_LP()
    @test isequal(cp, cp2)
    cp2.S[1] = 1
    @test !isequal(cp, cp2)
    @test isequal(cp, copy(cp))

    cp = test_coupledLP()
    @test nCouplingConstraints(cp) == 2000
end
