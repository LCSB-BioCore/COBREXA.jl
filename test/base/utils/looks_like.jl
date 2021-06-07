@testset "Find exchange reactions and metabolites CoreModel" begin
    cp = test_LP()
    @test isempty(find_exchange_reactions(cp))

    cp = test_simpleLP()
    @test isempty(find_exchange_reactions(cp))

    cp = CoreModel(
        [-1.0 -1 -2; 0 -1 0; 0 0 0],
        zeros(3),
        ones(3),
        ones(3),
        ones(3),
        ["EX_m1"; "r2"; "r3"],
        ["m1"; "m2"; "m3"],
    )
    @test find_exchange_reactions(cp) == ["EX_m1"]

    cp = CoreModel(
        [-1.0 0 0; 0 0 -1; 0 -1 0],
        zeros(3),
        ones(3),
        ones(3),
        ones(3),
        ["EX_m1"; "Exch_m3"; "Ex_m2"],
        ["m1"; "m2"; "m3"],
    )
    @test filter(looks_like_exchange_reaction, reactions(cp)) == ["EX_m1", "Exch_m3", "Ex_m2"]
    @test filter(x -> looks_like_exchange_reaction(x; ex_prefixes = ["Exch_"]), reactions(cp)) == ["Exch_m3"]

    # this is originally the "toyModel1.mat"
    cp = test_toyModel()

    @test find_exchange_reactions(cp) == ["EX_m1(e)", "EX_m3(e)", "EX_biomass(e)", "biomass1"]
    @test find_exchange_reactions(cp, exclude_biomass = true) == ["EX_m1(e)", "EX_m3(e)"]
end

@testset "Find exchange reactions and metabolites: all models except core" begin
    model = load_model(model_paths["e_coli_core.json"])
    @test length(filter(looks_like_exchange_reaction, reactions(model))) == 21
    
    model = load_model(model_paths["e_coli_core.xml"])
    @test length(filter(looks_like_exchange_reaction, reactions(model))) == 21
   
    model = load_model(model_paths["e_coli_core.mat"])
    @test length(filter(looks_like_exchange_reaction, reactions(model))) == 21
    
    model = convert(StandardModel, model)
    @test length(filter(looks_like_exchange_reaction, reactions(model))) == 21
   
    model = convert(CoreModelCoupled, model)
    @test length(filter(looks_like_exchange_reaction, reactions(model))) == 21
end
