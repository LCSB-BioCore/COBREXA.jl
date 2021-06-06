@testset "Change bounds" begin
    cp = test_LP()
    change_bounds!(cp, [3; 1], xl = [-10.0; -20], xu = [10.0; 20])
    @test cp.xl == [-20; 1; -10]
    @test cp.xu == [20; 1; 10]
    change_bounds!(
        cp,
        ["gibberish1"; "r3"; "r1"; "gibberish2"],
        xl = [0; -30.0; -40; 0],
        xu = [0; 30.0; 40; 0],
    )
    @test cp.xl == [-40; 1; -30]
    @test cp.xu == [40; 1; 30]
    change_bounds!(cp, ["r1"; "r3"], xl = [-50.0; -60])
    @test cp.xl == [-50; 1; -60]
end

@testset "Verify consistency" begin
    cp = test_LP()
    @test size(cp.S) == (4, 3)
    (new_reactions, new_mets) = verify_consistency(
        cp,
        reshape(cp.S[:, end], :, 1),
        [1.0, 2.0, 3.0, 4.0],
        [2.0],
        [-1.0],
        [1.0],
        ["r4"],
        ["m1", "m2", "m3", "m6"],
        [1],
        [4],
    )
    @test new_reactions == [1]
    @test new_mets == [4]

    (new_reactions, new_mets) = verify_consistency(
        cp,
        reshape(cp.S[:, end], :, 1),
        [1.0, 2.0, 3.0, 4.0],
        [2.0],
        [-1.0],
        [1.0],
        ["r1"],
        ["m1", "m2", "m3", "m6"],
        [],
        [4],
    )
    @test new_reactions == []
    @test new_mets == [4]
end

@testset "Add reactions (checking existence and consistency)" begin
    cp = test_LP()
    @test size(cp.S) == (4, 3)
    (new_cp, new_reactions, new_mets) = add_reactions(
        cp,
        cp.S[:, end],
        [1.0, 2.0, 3.0, 4.0],
        2.0,
        -1.0,
        1.0,
        check_consistency = true,
    )
    @test n_reactions(cp) + 1 == n_reactions(new_cp)

    (new_cp, new_reactions, new_mets) = add_reactions(
        cp,
        cp.S[:, end],
        [1.0, 2.0, 3.0, 4.0],
        2.0,
        -1.0,
        1.0,
        "r1",
        ["m1", "m2", "m3", "m6"],
        check_consistency = true,
    )
    @test n_reactions(cp) == n_reactions(new_cp)
    @test n_metabolites(cp) + 1 == n_metabolites(new_cp)
end

@testset "Add reactions" begin
    cp = test_LP()
    @test size(cp.S) == (4, 3)
    cp = add_reactions(cp, 2.0 * ones(4), 3 .* ones(4), 2.0, -1.0, 1.0)
    @test size(cp.S) == (8, 4)
    cp = add_reactions(cp, 2.0 * ones(4, 1), 3 .* ones(4), 2 .* ones(1), -ones(1), ones(1))
    @test size(cp.S) == (12, 5)
    cp = add_reactions(
        cp,
        2.0 * ones(4, 10),
        3 .* ones(4),
        2 .* ones(10),
        -ones(10),
        ones(10),
    )
    @test size(cp.S) == (16, 15)

    cp = test_sparseLP()
    @test size(cp.S) == (4000, 3000)
    cp = add_reactions(cp, 2.0 * sprand(4000, 0.5), 3 .* sprand(4000, 0.5), 2.0, -1.0, 1.0)
    @test size(cp.S) == (8000, 3001)
    cp = add_reactions(
        cp,
        2.0 * sprand(4000, 1, 0.5),
        3 .* sprand(4000, 0.5),
        2 .* sprand(1, 0.5),
        -sprand(1, 0.5),
        sprand(1, 0.5),
    )
    @test size(cp.S) == (12000, 3002)
    cp = add_reactions(
        cp,
        2.0 * sprand(4000, 1000, 0.5),
        3 .* sprand(4000, 0.5),
        2 .* sprand(1000, 0.5),
        -sprand(1000, 0.5),
        sprand(1000, 0.5),
    )
    @test size(cp.S) == (16000, 4002)

    cp = test_sparseLP()
    @test size(cp.S) == (4000, 3000)
    cp = add_reactions(cp, 2.0 * ones(4000), 3 .* ones(4000), 2.0, -1.0, 1.0)
    @test size(cp.S) == (8000, 3001)
    cp = add_reactions(
        cp,
        2.0 * ones(4000, 1),
        3 .* ones(4000),
        2 .* ones(1),
        -ones(1),
        ones(1),
    )
    @test size(cp.S) == (12000, 3002)
    cp = add_reactions(
        cp,
        2.0 * ones(4000, 1000),
        3 .* ones(4000),
        2 .* ones(1000),
        -ones(1000),
        ones(1000),
    )
    @test size(cp.S) == (16000, 4002)

    # proper subset of existing metabolites
    cp = test_LP()
    new_cp = add_reactions(cp, [-1.0], zeros(1), 1.0, 0.0, 1.0, "r4", ["m1"])
    @test n_reactions(cp) + 1 == n_reactions(new_cp)

    @test_throws DimensionMismatch add_reactions(
        cp,
        2.0 * ones(4000, 1),
        3 .* ones(4000),
        2 .* ones(2),
        -ones(1),
        ones(1),
    )
end

@testset "Remove reactions" begin
    cp = test_LP()
    @test size(cp.S) == (4, 3)
    cp = remove_reactions(cp, 2)
    @test size(cp.S) == (0, 2)
    cp = remove_reactions(cp, [2, 1])
    @test size(cp.S) == (0, 0)

    cp = test_LP()
    @test size(cp.S) == (4, 3)
    cp = remove_reactions(cp, "r0")
    @test size(cp.S) == (4, 3)
    cp = remove_reactions(cp, "r1")
    @test size(cp.S) == (0, 2)
    cp = remove_reactions(cp, ["r2"])
    @test size(cp.S) == (0, 1)

    lp = CoreModel(
        [1.0 1 1 0; 1 1 1 0; 1 1 1 0; 0 0 0 1],
        collect(1.:4),
        collect(1.:4),
        collect(1.:4),
        collect(1.:4),
        ["r1"; "r2"; "r3"; "r4"],
        ["m1"; "m2"; "m3"; "m4"],
    )

    modLp = remove_reactions(lp, [4; 1])
    @test stoichiometry(modLp) == stoichiometry(lp)[1:3, 2:3]
    @test balance(modLp) == balance(lp)[1:3]
    @test objective(modLp) == objective(lp)[2:3]
    @test bounds(modLp)[1] == bounds(lp)[1][2:3]
    @test bounds(modLp)[2] == bounds(lp)[2][2:3]
    @test reactions(modLp) == reactions(lp)[2:3]
    @test metabolites(modLp) == metabolites(lp)[1:3]
end

@testset "Remove metabolites" begin
    model = load_model(CoreModel, model_paths["e_coli_core.json"])

    m1 = remove_metabolites(model, ["glc__D_e", "for_c"])
    m2 = remove_metabolites(model, "glc__D_e")
    m3 = remove_metabolites(model, indexin(["glc__D_e", "for_c"], metabolites(model)))
    m4 = remove_metabolites(model, first(indexin(["glc__D_e"], metabolites(model))))

    @test size(stoichiometry(m1)) == (70, 94)
    @test size(stoichiometry(m2)) == (71, 94)
    @test size(stoichiometry(m3)) == (70, 94)
    @test size(stoichiometry(m4)) == (71, 94)
    @test any(["glc__D_e", "for_c"] .∉ Ref(metabolites(m1)))
    @test any(["glc__D_e"] .∉ Ref(metabolites(m2)))
    @test any(["glc__D_e", "for_c"] .∉ Ref(metabolites(m3)))
    @test any(["glc__D_e"] .∉ Ref(metabolites(m4)))
end
