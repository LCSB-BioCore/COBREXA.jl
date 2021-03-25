@testset "Coupling constraints" begin
    cp = convert(CoupledLinearModel, test_LP())
    @test size(cp.lm.S) == (4, 3)
    @test size(stoichiometry(convert(LinearModel, cp))) == (4, 3)
    newCp = addCouplingConstraints(cp, stoichiometry(cp)[end, :], -1.0, 1.0)
    @test nCouplingConstraints(cp) + 1 == nCouplingConstraints(newCp)

    newCp = addCouplingConstraints(cp, stoichiometry(cp)[1:2, :], [-1.0; -1.0], [1.0; 1.0])
    @test nCouplingConstraints(cp) + 2 == nCouplingConstraints(newCp)

    nC = nCouplingConstraints(cp)
    addCouplingConstraints!(cp, stoichiometry(cp)[end, :], -1.0, 1.0)
    @test nC + 1 == nCouplingConstraints(cp)
    addCouplingConstraints!(cp, stoichiometry(cp)[1:2, :], [-1.0; -1.0], [1.0; 1.0])
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
    @test size(coupling(cp)) == (nC, nReactions(cp))
    @test nC - 1 == nCouplingConstraints(newCp)
    @test nCouplingConstraints(cp) == nC
    newCp = removeCouplingConstraints(cp, [1, 2])
    @test nC - 2 == nCouplingConstraints(newCp)
    newCp = removeCouplingConstraints(cp, Array(1:nCouplingConstraints(cp)))
    @test nCouplingConstraints(newCp) == 0
    @test nCouplingConstraints(cp) == nC

    cp = test_coupledLP()
    changeCouplingBounds!(cp, [3, 1], cl = [-10.0, -20], cu = [10.0, 20])
    cl, cu = couplingBounds(cp)
    @test cl[[1, 3]] == [-20, -10]
    @test cu[[1, 3]] == [20, 10]
    changeCouplingBounds!(cp, [1000, 1001], cl = [-50.0, -60.0])
    cl, cu = couplingBounds(cp)
    @test cl[[1000, 1001]] == [-50.0, -60.0]
end
