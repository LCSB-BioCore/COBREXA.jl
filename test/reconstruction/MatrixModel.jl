@testset "Change bounds" begin
    cp = test_LP()

    change_bound!(cp, 1, lower = -10, upper = 10)
    @test cp.xl[1] == -10
    @test cp.xu[1] == 10
    change_bounds!(cp, [1, 2]; lower = [-11, -12.2], upper = [11, 23.0])
    @test cp.xl[2] == -12.2
    @test cp.xu[1] == 11

    change_bound!(cp, "r1", lower = -101, upper = 101)
    @test cp.xl[1] == -101
    @test cp.xu[1] == 101
    change_bounds!(cp, ["r1", "r2"]; lower = [-113, -12.23], upper = [114, 233.0])
    @test cp.xl[2] == -12.23
    @test cp.xu[1] == 114
    change_bounds!(cp, ["r1", "r2"]; lower = [-114, nothing], upper = [nothing, 2333.0])
    @test cp.xl[1] == -114
    @test cp.xl[2] == -12.23
    @test cp.xu[1] == 114
    @test cp.xu[2] == 2333

    new_model = change_bound(cp, 1, lower = -10, upper = 10)
    @test new_model.xl[1] == -10
    @test new_model.xu[1] == 10
    new_model = change_bounds(cp, [1, 2]; lower = [-11, -12.2], upper = [11, 23.0])
    @test new_model.xl[2] == -12.2
    @test new_model.xu[1] == 11

    new_model = change_bound(cp, "r1", lower = -101, upper = 101)
    @test new_model.xl[1] == -101
    @test new_model.xu[1] == 101
    new_model = change_bounds(cp, ["r1", "r2"]; lower = [-113, -12.23], upper = [113, 1000])
    @test new_model.xl[2] == -12.23
    @test new_model.xu[1] == 113
    new_model =
        change_bounds(cp, ["r1", "r2"]; lower = [nothing, -10], upper = [110, nothing])
    @test new_model.xl[1] == -114
    @test new_model.xl[2] == -10
    @test new_model.xu[1] == 110
    @test new_model.xu[2] == 2333
end

@testset "Verify consistency" begin
    cp = test_LP()
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
    (new_cp, new_reactions, new_mets) = add_reactions(
        cp,
        cp.S[:, end],
        [1.0, 2.0, 3.0, 4.0],
        2.0,
        -1.0,
        1.0,
        check_consistency = true,
    )
    @test n_variables(cp) + 1 == n_variables(new_cp)

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
    @test n_variables(cp) == n_variables(new_cp)
    @test n_metabolites(cp) + 1 == n_metabolites(new_cp)
end

@testset "Add reactions" begin
    cp = test_LP()
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
    @test n_variables(cp) + 1 == n_variables(new_cp)

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
    cp = remove_reaction(cp, 2)
    @test size(cp.S) == (4, 2)
    cp = remove_reactions(cp, [2, 1])
    @test size(cp.S) == (4, 0)

    cp = test_LP()
    cp = remove_reaction(cp, "r1")
    @test size(cp.S) == (4, 2)
    cp = remove_reactions(cp, ["r2"])
    @test size(cp.S) == (4, 1)

    lp = MatrixModel(
        [1.0 1 1 0; 1 1 1 0; 1 1 1 0; 0 0 0 1],
        collect(1.0:4),
        collect(1.0:4),
        collect(1.0:4),
        collect(1.0:4),
        ["r1"; "r2"; "r3"; "r4"],
        ["m1"; "m2"; "m3"; "m4"],
    )

    modLp = remove_reactions(lp, [4; 1])
    @test stoichiometry(modLp) == stoichiometry(lp)[:, 2:3]
    @test balance(modLp) == balance(lp)
    @test objective(modLp) == objective(lp)[2:3]
    @test bounds(modLp)[1] == bounds(lp)[1][2:3]
    @test bounds(modLp)[2] == bounds(lp)[2][2:3]
    @test reactions(modLp) == reactions(lp)[2:3]
    @test metabolites(modLp) == metabolites(lp)
end

@testset "Remove metabolites" begin
    model = load_model(MatrixModel, model_paths["e_coli_core.json"])

    m1 = remove_metabolites(model, ["glc__D_e", "for_c"])
    m2 = remove_metabolite(model, "glc__D_e")
    m3 = remove_metabolites(model, Int.(indexin(["glc__D_e", "for_c"], metabolites(model))))
    m4 = remove_metabolite(model, first(indexin(["glc__D_e"], metabolites(model))))

    @test size(stoichiometry(m1)) == (70, 90)
    @test size(stoichiometry(m2)) == (71, 93)
    @test size(stoichiometry(m3)) == (70, 90)
    @test size(stoichiometry(m4)) == (71, 93)
    @test all((!in(metabolites(m1))).(["glc__D_e", "for_c"]))
    @test !(["glc__D_e"] in metabolites(m2))
    @test all((!in(metabolites(m3))).(["glc__D_e", "for_c"]))
    @test !(["glc__D_e"] in metabolites(m4))
end

@testset "Core in place modifications" begin
    toymodel = test_toyModel()

    rxn1 = Reaction("nr1"; metabolites = Dict("m1[c]" => -1, "m3[c]" => 1))
    rxn2 = Reaction("nr2"; metabolites = Dict("m1[c]" => -1, "m2[c]" => 1))
    rxn3 = Reaction("nr3"; metabolites = Dict("m2[c]" => -1, "m3[c]" => 1))
    rxn3.lower_bound = 10

    add_reaction!(toymodel, rxn1)
    @test toymodel.S[1, 8] == -1
    @test toymodel.S[2, 8] == 1
    @test all(toymodel.S[3:end, 8] .== 0)
    @test size(toymodel.S) == (6, 8)
    @test toymodel.rxns[end] == "nr1"

    add_reactions!(toymodel, [rxn2, rxn3])
    @test size(toymodel.S) == (6, 10)
    @test toymodel.xl[end] == 10

    change_objective!(toymodel, "nr1")
    @test objective(toymodel)[8] == 1.0
    @test objective(toymodel)[7] == 0.0
end
