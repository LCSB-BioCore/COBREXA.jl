@testset "Find exchange reactions and metabolites" begin
    cp = test_LP()
    @test isempty(find_exchange_reactions(cp))
    @test isempty(find_exchange_metabolites(cp))

    cp = test_simpleLP()
    @test isempty(find_exchange_reactions(cp))
    @test isempty(find_exchange_metabolites(cp))

    cp = CoreModel(
        [-1.0 -1 -2; 0 -1 0; 0 0 0],
        zeros(3),
        ones(3),
        ones(3),
        ones(3),
        ["EX_m1"; "r2"; "r3"],
        ["m1"; "m2"; "m3"],
    )
    @test find_exchange_reactions(cp) == [1]

    cp = CoreModel(
        [-1.0 0 0; 0 0 -1; 0 -1 0],
        zeros(3),
        ones(3),
        ones(3),
        ones(3),
        ["EX_m1"; "Exch_m3"; "Ex_m2"],
        ["m1"; "m2"; "m3"],
    )
    @test find_exchange_reactions(cp) == [1; 2; 3]
    @test find_exchange_metabolites(cp) == [1; 3; 2]
    @test find_exchange_reactions(cp, exc_prefs = ["Exch_"]) == [2]
    @test find_exchange_metabolites(cp, exc_prefs = ["Exch_"]) == [3]

    # this is originally the "toyModel1.mat"
    cp = test_toyModel()

    @test find_exchange_reactions(cp) == [4; 5; 6]
    @test find_exchange_metabolites(cp) == [4; 5; 6]
    @test find_exchange_reactions(cp, exclude_biomass = true) == [4; 5]
    @test find_exchange_metabolites(cp, exclude_biomass = true) == [4; 5]
    @test find_exchange_reactions(cp, exclude_biomass = true, biomass_str = "biom") ==
          [4; 5]
    @test find_exchange_metabolites(cp, exclude_biomass = true, biomass_str = "biom") ==
          [4; 5]
end


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
end
