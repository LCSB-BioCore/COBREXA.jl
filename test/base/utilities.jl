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

@testset "Allocating reactions" begin
    @test COBREXA.allocateReacs([1], 1) == [[1]]
    @test COBREXA.allocateReacs(collect(1:3), 2) == [[1;2], [3]]
    @test COBREXA.allocateReacs(collect(1:8), 4) == [[1;2], [3;4], [5;6], [7;8]]
    @test COBREXA.allocateReacs(collect(1:11), 3) == [[1;2;3;4], [5;6;7;8], [9;10;11]]
end
