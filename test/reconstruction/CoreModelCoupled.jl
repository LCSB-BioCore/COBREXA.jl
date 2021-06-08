@testset "Coupling constraints" begin
    cp = convert(CoreModelCoupled, test_LP())
    @test size(cp.lm.S) == (4, 3)
    @test size(stoichiometry(convert(CoreModel, cp))) == (4, 3)
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

@testset "Add reactions" begin
    cp = convert(CoreModelCoupled, test_LP())
    cp = add_coupling_constraints(cp, stoichiometry(cp)[end, :], -1.0, 1.0)

    new_cp = add_reactions(cp, 2.0 * ones(4), 3 .* ones(4), 2.0, -1.0, 1.0)
    @test new_cp isa CoreModelCoupled
    @test cp.C == new_cp.C[:, 1:end-1]
    @test cp.cl == new_cp.cl
    @test cp.cu == new_cp.cu

    new_cp = add_reactions(cp, 2.0 * ones(4), 3 .* ones(4), 2.0, -1.0, 1.0, "r4", ["m$i" for i in 1:4])
    @test cp.C == new_cp.C[:, 1:end-1]
    @test cp.cl == new_cp.cl
    @test cp.cu == new_cp.cu

    new_cp = add_reactions(
        cp,
        2.0 * ones(4, 10),
        3 .* ones(4),
        2 .* ones(10),
        -ones(10),
        ones(10),
    )
    @test cp.C == new_cp.C[:, 1:end-10]
    @test cp.cl == new_cp.cl
    @test cp.cu == new_cp.cu

    new_cp = add_reactions(
        cp,
        2.0 * ones(4, 10),
        3 .* ones(4),
        2 .* ones(10),
        -ones(10),
        ones(10),
        ["r$i" for i in 1:10],
        ["m$i" for i in 1:4]
    )
    @test cp.C == new_cp.C[:, 1:end-7] # 3 reactions were already present
    @test cp.cl == new_cp.cl
    @test cp.cu == new_cp.cu

    new_cp = add_reactions(
            cp,
            2.0 * sprand(4000, 0.5),
            3 .* sprand(4000, 0.5),
            2.0,
            -1.0,
            1.0
        )
    @test cp.C == new_cp.C[:, 1:end-1]
    @test cp.cl == new_cp.cl
    @test cp.cu == new_cp.cu

    cm = CoreModel(
            2.0 * ones(4, 10),
            3 .* ones(4),
            2 .* ones(10),
            -ones(10),
            ones(10),
            ["r$i" for i in 1:10],
            ["m$i" for i in 1:4]
        )
    new_cp = add_reactions(cp, cm)
    @test cp.C == new_cp.C[:, 1:end-7] # 3 reactions were already present
    @test cp.cl == new_cp.cl
    @test cp.cu == new_cp.cu
end

@testset "Remove reactions" begin
    cp = convert(CoreModelCoupled, test_LP())
    cp = add_coupling_constraints(cp, 1. .* collect(1:n_reactions(cp)), -1.0, 1.0)

    new_cp = remove_reactions(cp, [3;2])
    @test new_cp isa CoreModelCoupled
    @test new_cp.C[:] == cp.C[:, 1] # because cp.C[:, 1] comes out as a Vector
    @test new_cp.cl == cp.cl
    @test new_cp.cu == cp.cu

    new_cp = remove_reactions(cp, 2)
    @test new_cp.C == cp.C[:, [1; 3]]
    @test new_cp.cl == cp.cl
    @test new_cp.cu == cp.cu

    new_cp = remove_reactions(cp, "r1")
    @test new_cp.C == cp.C[:, 2:3]
    @test new_cp.cl == cp.cl
    @test new_cp.cu == cp.cu

    new_cp = remove_reactions(cp, ["r1"; "r3"])
    @test new_cp.C[:] == cp.C[:, 2]
    @test new_cp.cl == cp.cl
    @test new_cp.cu == cp.cu

    new_cp = remove_reactions(cp, [1;4])
    @test new_cp.C == cp.C[:, 2:3]

    new_cp = remove_reactions(cp, "r4")
    @test new_cp.C == cp.C

    new_cp = remove_reactions(cp, [1;1;2])
    @test new_cp.C[:] == cp.C[:, 3]
end

@testset "Change bounds" begin
    cp = convert(CoreModelCoupled, test_LP())
    change_bounds!(cp, [3; 1], xl = [-10.0; -20], xu = [10.0; 20])
    @test cp isa CoreModelCoupled
    @test cp.lm.xl == [-20; 1; -10]
    @test cp.lm.xu == [20; 1; 10]
    change_bounds!(
        cp,
        ["gibberish1"; "r3"; "r1"; "gibberish2"],
        xl = [0; -30.0; -40; 0],
        xu = [0; 30.0; 40; 0],
    )
    @test cp.lm.xl == [-40; 1; -30]
    @test cp.lm.xu == [40; 1; 30]
    change_bounds!(cp, ["r1"; "r3"], xl = [-50.0; -60])
    @test cp.lm.xl == [-50; 1; -60]
end
