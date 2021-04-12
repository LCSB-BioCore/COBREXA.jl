@testset "Coupling constraints" begin
    cp = convert(CoupledLinearModel, test_LP())
    @test size(cp.lm.S) == (4, 3)
    @test size(stoichiometry(convert(LinearModel, cp))) == (4, 3)
    new_cp = add_coupling_constraints(cp, stoichiometry(cp)[end, :], -1.0, 1.0)
    @test n_coupling_constraints(cp) + 1 == n_coupling_constraints(new_cp)

    new_cp =
        add_coupling_constraints(cp, stoichiometry(cp)[1:2, :], [-1.0; -1.0], [1.0; 1.0])
    @test n_coupling_constraints(cp) + 2 == n_coupling_constraints(new_cp)

    n_c = n_coupling_constraints(cp)
    add_coupling_constraints!(cp, stoichiometry(cp)[end, :], -1.0, 1.0)
    @test n_c + 1 == n_coupling_constraints(cp)
    add_coupling_constraints!(cp, stoichiometry(cp)[1:2, :], [-1.0; -1.0], [1.0; 1.0])
    @test n_c + 3 == n_coupling_constraints(cp)

    n_c = n_coupling_constraints(cp)
    remove_coupling_constraints!(cp, 1)
    @test n_c - 1 == n_coupling_constraints(cp)
    remove_coupling_constraints!(cp, [1, 2])
    @test n_c - 3 == n_coupling_constraints(cp)
    @test n_coupling_constraints(cp) == 0

    cp = test_coupledLP()
    n_c = n_coupling_constraints(cp)
    new_cp = remove_coupling_constraints(cp, 1)
    @test size(coupling(cp)) == (n_c, n_reactions(cp))
    @test n_c - 1 == n_coupling_constraints(new_cp)
    @test n_coupling_constraints(cp) == n_c
    new_cp = remove_coupling_constraints(cp, [1, 2])
    @test n_c - 2 == n_coupling_constraints(new_cp)
    new_cp = remove_coupling_constraints(cp, Vector(1:n_coupling_constraints(cp)))
    @test n_coupling_constraints(new_cp) == 0
    @test n_coupling_constraints(cp) == n_c

    cp = test_coupledLP()
    change_coupling_bounds!(cp, [3, 1], cl = [-10.0, -20], cu = [10.0, 20])
    cl, cu = coupling_bounds(cp)
    @test cl[[1, 3]] == [-20, -10]
    @test cu[[1, 3]] == [20, 10]
    change_coupling_bounds!(cp, [1000, 1001], cl = [-50.0, -60.0])
    cl, cu = coupling_bounds(cp)
    @test cl[[1000, 1001]] == [-50.0, -60.0]
end
