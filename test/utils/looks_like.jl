@testset "Looks like in MatrixModel, detailed test" begin
    cp = test_LP()
    @test isempty(filter(looks_like_exchange_reaction, variables(cp)))

    cp = test_simpleLP()
    @test isempty(filter(looks_like_exchange_reaction, variables(cp)))

    cp = MatrixModel(
        [-1.0 -1 -2; 0 -1 0; 0 0 0],
        zeros(3),
        ones(3),
        ones(3),
        ones(3),
        ["EX_m1"; "r2"; "r3"],
        ["m1"; "m2"; "m3"],
    )
    @test filter(looks_like_exchange_reaction, variables(cp)) == ["EX_m1"]

    cp = MatrixModel(
        [-1.0 0 0; 0 0 -1; 0 -1 0],
        zeros(3),
        ones(3),
        ones(3),
        ones(3),
        ["EX_m1"; "Exch_m3"; "Ex_m2"],
        ["m1"; "m2"; "m3_e"],
    )
    @test filter(looks_like_exchange_reaction, variables(cp)) ==
          ["EX_m1", "Exch_m3", "Ex_m2"]
    @test filter(
        x -> looks_like_exchange_reaction(x; exchange_prefixes = ["Exch_"]),
        variables(cp),
    ) == ["Exch_m3"]

    # this is originally the "toyModel1.mat"
    cp = test_toyModel()

    @test filter(looks_like_exchange_reaction, variables(cp)) ==
          ["EX_m1(e)", "EX_m3(e)", "EX_biomass(e)"]
    @test filter(
        x -> looks_like_exchange_reaction(x; exclude_biomass = true),
        variables(cp),
    ) == ["EX_m1(e)", "EX_m3(e)"]
    @test filter(looks_like_extracellular_metabolite, metabolites(cp)) == ["m1[e]", "m3[e]"]
    @test filter(looks_like_biomass_reaction, variables(cp)) ==
          ["EX_biomass(e)", "biomass1"]
    @test filter(
        x -> looks_like_biomass_reaction(x; exclude_exchanges = true),
        variables(cp),
    ) == ["biomass1"]
end

@testset "Looks like functions, basic" begin
    model = load_model(model_paths["e_coli_core.json"])
    @test length(filter(looks_like_exchange_reaction, variables(model))) == 20
    @test length(filter(looks_like_biomass_reaction, variables(model))) == 1
    @test length(filter(looks_like_extracellular_metabolite, metabolites(model))) == 20

    model = load_model(model_paths["e_coli_core.xml"])
    @test length(filter(looks_like_exchange_reaction, variables(model))) == 20
    @test length(filter(looks_like_biomass_reaction, variables(model))) == 1
    @test length(filter(looks_like_extracellular_metabolite, metabolites(model))) == 20

    model = load_model(model_paths["e_coli_core.mat"])
    @test length(filter(looks_like_exchange_reaction, variables(model))) == 20
    @test length(filter(looks_like_biomass_reaction, variables(model))) == 1
    @test length(filter(looks_like_extracellular_metabolite, metabolites(model))) == 20

    model = convert(ObjectModel, model)
    @test length(filter(looks_like_exchange_reaction, variables(model))) == 20
    @test length(filter(looks_like_biomass_reaction, variables(model))) == 1
    @test length(filter(looks_like_extracellular_metabolite, metabolites(model))) == 20

    model = convert(MatrixModelWithCoupling, model)
    @test length(filter(looks_like_exchange_reaction, variables(model))) == 20
    @test length(filter(looks_like_biomass_reaction, variables(model))) == 1
    @test length(filter(looks_like_extracellular_metabolite, metabolites(model))) == 20
end

@testset "Ontology usage in is_xxx_reaction" begin
    model = load_model(ObjectModel, model_paths["e_coli_core.json"])

    # macro generated, so only test positive and negative case
    @test !is_biomass_reaction(model, "PFL")
    @test is_biomass_reaction(model, "BIOMASS_Ecoli_core_w_GAM")

    @test isa_reaction(model, "PFL")
    @test_throws KeyError isa_reaction(model, "atp_c") 
    @test isa_gene(model, "b2464")
    @test_throws KeyError isa_gene(model, "atp_c")
    @test isa_metabolite(model, "atp_c")
    @test_throws KeyError isa_metabolite(model, "PFL")
end
