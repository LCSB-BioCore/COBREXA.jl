@testset "Remove reactions" begin
lp = LinearModel([1. 1 1 0; 1 1 1 0; 1 1 1 0; 0 0 0 1],
                 zeros(4),
                 zeros(4),
                 zeros(4),
                 zeros(4),
                 ["r1"; "r2"; "r3"; "r4"],
                 ["m1"; "m2"; "m3"; "m4"])
end

@testset "Find exchange reactions and metabolites" begin
    cp = test_LP()
    @test isempty(findExchangeReactions(cp))
    @test isempty(findExchangeMetabolites(cp))

    cp = test_simpleLP()
    @test isempty(findExchangeReactions(cp))
    @test isempty(findExchangeMetabolites(cp))

    cp = LinearModel([-1. -1 -2; 0 -1 0; 0 0 0],
                        zeros(3),
                        ones(3),
                        ones(3),
                        ones(3),
                        ["EX_m1"; "r2"; "r3"],
                        ["m1"; "m2"; "m3"])
    @test findExchangeReactions(cp)==[1]

    cp = LinearModel([-1. 0 0; 0 0 -1; 0 -1 0],
                        zeros(3),
                        ones(3),
                        ones(3),
                        ones(3),
                        ["EX_m1"; "Exch_m3"; "Ex_m2"],
                        ["m1"; "m2"; "m3"])
    @test findExchangeReactions(cp)==[1;2;3]
    @test findExchangeMetabolites(cp)==[1;3;2]
    @test findExchangeReactions(cp, excPrefs=["Exch_"])==[2]
    @test findExchangeMetabolites(cp, excPrefs=["Exch_"])==[3]

    cp = loadModel(joinpath("data", "toyModel1.mat"), "model")
    @test findExchangeReactions(cp)==[4;5;6]
    @test findExchangeMetabolites(cp)==[4;5;6]
    @test findExchangeReactions(cp, excludeBiomass=true)==[4;5]
    @test findExchangeMetabolites(cp, excludeBiomass=true)==[4;5]
    @test findExchangeReactions(cp, excludeBiomass=true, biomassStr="biom")==[4;5]
    @test findExchangeMetabolites(cp, excludeBiomass=true, biomassStr="biom")==[4;5]
end


@testset "Change bounds" begin
  cp = test_LP()
  changeBounds!(cp, [3; 1], xl=[-10.; -20], xu=[10.; 20])
  @test cp.xl == [-20; 1; -10]
  @test cp.xu == [20; 1; 10]
  changeBounds!(cp, ["gibberish1"; "r3"; "r1"; "gibberish2"], xl=[0; -30.; -40; 0], xu=[0; 30.; 40; 0])
  @test cp.xl == [-40; 1; -30]
  @test cp.xu == [40; 1; 30]
  changeBounds!(cp, ["r1"; "r3"], xl=[-50.; -60])
  @test cp.xl == [-50; 1; -60]
end

@testset "Verify consistency" begin
    cp = test_LP()
    @test size(cp.S) == (4, 3)
    (newReactions, newMets) = verifyConsistency(cp, reshape(cp.S[:, end], :, 1),
                                                [1., 2., 3., 4.], [2.], [-1.], [1.],
                                                ["r4"], ["m1", "m2", "m3", "m6"],
                                                [1], [4])
    @test newReactions == [1]
    @test newMets == [4]

    (newReactions, newMets) = verifyConsistency(cp, reshape(cp.S[:, end], :, 1),
                                                [1., 2., 3., 4.], [2.], [-1.], [1.],
                                                ["r1"], ["m1", "m2", "m3", "m6"],
                                                [], [4])
    @test newReactions == []
    @test newMets == [4]
end



@testset "Add reactions (checking existence and consistency)" begin
    cp = test_LP()
    @test size(cp.S) == (4, 3)
    (newCp, newReactions, newMets) = addReactions(cp, cp.S[:, end], [1., 2., 3., 4.], 2., -1., 1., checkConsistency=true)
    @test nReactions(cp) + 1 == nReactions(newCp)

    (newCp, newReactions, newMets) = addReactions(cp, cp.S[:, end], [1., 2., 3., 4.], 2., -1., 1., "r1", ["m1", "m2", "m3", "m6"], checkConsistency=true)
    @test nReactions(cp) == nReactions(newCp)
    @test nMetabolites(cp) + 1 == nMetabolites(newCp)
end

@testset "Add reactions" begin
    cp = test_LP()
    @test size(cp.S) == (4, 3)
    cp = addReactions(cp, 2. * ones(4), 3 .* ones(4), 2., -1., 1.)
    @test size(cp.S) == (8, 4)
    cp = addReactions(cp, 2. * ones(4, 1), 3 .* ones(4), 2 .* ones(1), -ones(1), ones(1))
    @test size(cp.S) == (12, 5)
    cp = addReactions(cp, 2. * ones(4, 10), 3 .* ones(4), 2 .* ones(10), -ones(10), ones(10))
    @test size(cp.S) == (16, 15)

    cp = test_sparseLP()
    @test size(cp.S) == (4000, 3000)
    cp = addReactions(cp, 2. * sprand(4000, 0.5), 3 .* sprand(4000, 0.5), 2., -1., 1.)
    @test size(cp.S) == (8000, 3001)
    cp = addReactions(cp, 2. * sprand(4000, 1, 0.5), 3 .* sprand(4000, 0.5), 2 .* sprand(1, 0.5), -sprand(1, 0.5), sprand(1, 0.5))
    @test size(cp.S) == (12000, 3002)
    cp = addReactions(cp, 2. * sprand(4000, 1000, 0.5), 3 .* sprand(4000, 0.5), 2 .* sprand(1000, 0.5), -sprand(1000, 0.5), sprand(1000, 0.5))
    @test size(cp.S) == (16000, 4002)

    cp = test_sparseLP()
    @test size(cp.S) == (4000, 3000)
    cp = addReactions(cp, 2. * ones(4000), 3 .* ones(4000), 2., -1., 1.)
    @test size(cp.S) == (8000, 3001)
    cp = addReactions(cp, 2. * ones(4000, 1), 3 .* ones(4000), 2 .* ones(1), -ones(1), ones(1))
    @test size(cp.S) == (12000, 3002)
    cp = addReactions(cp, 2. * ones(4000, 1000), 3 .* ones(4000), 2 .* ones(1000), -ones(1000), ones(1000))
    @test size(cp.S) == (16000, 4002)

    # proper subset of existing metabolites
    cp = test_LP()
    newCp = addReactions(cp, [-1.], zeros(1), 1., 0., 1., "r4", ["m1"])
    @test nReactions(cp) + 1 == nReactions(newCp)

    @test_throws DimensionMismatch addReactions(cp, 2. * ones(4000, 1), 3 .* ones(4000), 2 .* ones(2), -ones(1), ones(1))
end


@testset "Remove reactions" begin
    cp = test_LP()
    @test size(cp.S) == (4, 3)
    cp = removeReactions(cp, 2)
    @test size(cp.S) == (0, 2)
    cp = removeReactions(cp, [2, 1])
    @test size(cp.S) == (0, 0)

    cp = test_LP()
    @test size(cp.S) == (4, 3)
    cp = removeReactions(cp, "r0")
    @test size(cp.S) == (4, 3)
    cp = removeReactions(cp, "r1")
    @test size(cp.S) == (0, 2)
    cp = removeReactions(cp, ["r2"])
    @test size(cp.S) == (0, 1)
end
