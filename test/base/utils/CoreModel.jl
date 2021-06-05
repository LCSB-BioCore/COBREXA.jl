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

@testset "CoreModel utilities" begin
    cp = test_LP()
    @test n_reactions(cp) == 3
    @test n_metabolites(cp) == 4
    @test n_coupling_constraints(cp) == 0

    cp2 = test_LP()
    @test isequal(cp, cp2)
    cp2.S[1] = 1
    @test !isequal(cp, cp2)
    @test isequal(cp, copy(cp))

    cp = test_coupledLP()
    @test n_coupling_constraints(cp) == 2000
    @test isequal(cp, copy(cp))
end
