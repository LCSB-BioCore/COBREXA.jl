@testset "Looks like in CoreModel, detailed test" begin
    cp = test_LP()
    @test isempty(filter(looks_like_exchange_reaction, reactions(cp)))

    cp = test_simpleLP()
    @test isempty(filter(looks_like_exchange_reaction, reactions(cp)))

    cp = CoreModel(
        [-1.0 -1 -2; 0 -1 0; 0 0 0],
        zeros(3),
        ones(3),
        ones(3),
        ones(3),
        ["EX_m1"; "r2"; "r3"],
        ["m1"; "m2"; "m3"],
    )
    @test filter(looks_like_exchange_reaction, reactions(cp)) == ["EX_m1"]

    cp = CoreModel(
        [-1.0 0 0; 0 0 -1; 0 -1 0],
        zeros(3),
        ones(3),
        ones(3),
        ones(3),
        ["EX_m1"; "Exch_m3"; "Ex_m2"],
        ["m1"; "m2"; "m3_e"],
    )
    @test filter(looks_like_exchange_reaction, reactions(cp)) ==
          ["EX_m1", "Exch_m3", "Ex_m2"]
    @test filter(
        x -> looks_like_exchange_reaction(x; exchange_prefixes = ["Exch_"]),
        reactions(cp),
    ) == ["Exch_m3"]

    # this is originally the "toyModel1.mat"
    cp = test_toyModel()

    @test filter(looks_like_exchange_reaction, reactions(cp)) ==
          ["EX_m1(e)", "EX_m3(e)", "EX_biomass(e)"]
    @test filter(
        x -> looks_like_exchange_reaction(x; exclude_biomass = true),
        reactions(cp),
    ) == ["EX_m1(e)", "EX_m3(e)"]
    @test filter(looks_like_extracellular_metabolite, metabolites(cp)) == ["m1[e]", "m3[e]"]
    @test filter(looks_like_biomass_reaction, reactions(cp)) ==
          ["EX_biomass(e)", "biomass1"]
    @test filter(
        x -> looks_like_biomass_reaction(x; exclude_exchanges = true),
        reactions(cp),
    ) == ["biomass1"]
end

@testset "Looks like functions, basic" begin
    model = load_model(model_paths["e_coli_core.json"])
    @test length(filter(looks_like_exchange_reaction, reactions(model))) == 20
    @test length(filter(looks_like_biomass_reaction, reactions(model))) == 1
    @test length(filter(looks_like_extracellular_metabolite, metabolites(model))) == 20

    model = load_model(model_paths["e_coli_core.xml"])
    @test length(filter(looks_like_exchange_reaction, reactions(model))) == 20
    @test length(filter(looks_like_biomass_reaction, reactions(model))) == 1
    @test length(filter(looks_like_extracellular_metabolite, metabolites(model))) == 20

    model = load_model(model_paths["e_coli_core.mat"])
    @test length(filter(looks_like_exchange_reaction, reactions(model))) == 20
    @test length(filter(looks_like_biomass_reaction, reactions(model))) == 1
    @test length(filter(looks_like_extracellular_metabolite, metabolites(model))) == 20

    model = convert(StandardModel, model)
    @test length(filter(looks_like_exchange_reaction, reactions(model))) == 20
    @test length(filter(looks_like_biomass_reaction, reactions(model))) == 1
    @test length(filter(looks_like_extracellular_metabolite, metabolites(model))) == 20

    model = convert(CoreModelCoupled, model)
    @test length(filter(looks_like_exchange_reaction, reactions(model))) == 20
    @test length(filter(looks_like_biomass_reaction, reactions(model))) == 1
    @test length(filter(looks_like_extracellular_metabolite, metabolites(model))) == 20
end

@testset "Ontology usage in is_xxx_reaction" begin 
    model = load_model(StandardModel, model_paths["e_coli_core.json"])

    # macro generated, so only test positive and negative case
    @test !is_biomass_reaction(model, "PFL")
    @test is_biomass_reaction(model, "BIOMASS_Ecoli_core_w_GAM")
end