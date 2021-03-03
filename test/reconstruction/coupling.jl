@testset "Coupling constraints" begin
    cp = test_LP()
    @test size(cp.S) == (4, 3)
    newCp = addCouplingConstraints(cp, cp.S[end, :], -1.0, 1.0)
    @test nCouplingConstraints(cp) + 1 == nCouplingConstraints(newCp)

    newCp = addCouplingConstraints(cp, cp.S[1:2, :], [-1.0; -1.0], [1.0; 1.0])
    @test nCouplingConstraints(cp) + 2 == nCouplingConstraints(newCp)

    nC = nCouplingConstraints(cp)
    addCouplingConstraints!(cp, cp.S[end, :], -1.0, 1.0)
    @test nC + 1 == nCouplingConstraints(cp)
    addCouplingConstraints!(cp, cp.S[1:2, :], [-1.0; -1.0], [1.0; 1.0])
    @test nC + 3 == nCouplingConstraints(cp)

    nC = nCouplingConstraints(cp)
    removeCouplingConstraints!(cp, 1)
    @test nC - 1 == nCouplingConstraints(cp)
    removeCouplingConstraints!(cp, [1, 2])
    @test nC - 3 == nCouplingConstraints(cp)
    @test nCouplingConstraints(cp) == 0

    cp = test_coupledLP()
    nC = nCouplingConstraints(cp)
    newCp = removeCouplingConstraints(cp, 1)
    @test nC - 1 == nCouplingConstraints(newCp)
    @test nCouplingConstraints(cp) == nC
    newCp = removeCouplingConstraints(cp, [1, 2])
    @test nC - 2 == nCouplingConstraints(newCp)
    newCp = removeCouplingConstraints(cp, Array(1:nCouplingConstraints(cp)))
    @test nCouplingConstraints(newCp) == 0
    @test nCouplingConstraints(cp) == nC

    cp = test_coupledLP()
    changeCouplingBounds!(cp, [3, 1], cl = [-10.0, -20], cu = [10.0, 20])
    @test cp.cl[[1, 3]] == [-20, -10]
    @test cp.cu[[1, 3]] == [20, 10]
    changeCouplingBounds!(cp, [1000, 1001], cl = [-50.0, -60.0])
    @test cp.cl[[1000, 1001]] == [-50.0, -60.0]
end
