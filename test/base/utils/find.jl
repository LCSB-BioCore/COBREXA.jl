@testset "Find exchange reactions and metabolites CoreModel" begin
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
    @test find_exchange_reactions(cp) == ["EX_m1", "Exch_m3", "Ex_m2"]
    @test find_exchange_metabolites(cp)["Exch_m3"]["m3"] == -1.0 
    @test find_exchange_reactions(cp, ex_prefixes = ["Exch_"]) == ["Exch_m3"]
    @test find_exchange_metabolites(cp, ex_prefixes = ["Exch_"])["Exch_m3"]["m3"] == -1.0

    # this is originally the "toyModel1.mat"
    cp = test_toyModel()

    @test find_exchange_reactions(cp) == ["EX_m1(e)", "EX_m3(e)", "EX_biomass(e)", "biomass1"]
    @test find_exchange_metabolites(cp)["EX_biomass(e)"]["biomass[c]"] == -1.0
    @test find_exchange_reactions(cp, exclude_biomass = true) == ["EX_m1(e)", "EX_m3(e)"]
    @test find_exchange_metabolites(cp, exclude_biomass = true)["EX_m3(e)"]["m3[e]"] == -1.0
end

@testset "Find exchange reactions and metabolites: all models except core" begin
    model = load_model(model_paths["e_coli_core.json"])
    @test length(find_exchange_reactions(model)) == 21
    @test length(find_exchange_metabolites(model)["EX_fum_e"]) == 1
    
    model = load_model(model_paths["e_coli_core.xml"])
    @test length(find_exchange_reactions(model)) == 21
    @test length(find_exchange_metabolites(model)["R_EX_fum_e"]) == 1
   
    model = load_model(model_paths["e_coli_core.mat"])
    @test length(find_exchange_reactions(model)) == 21
    @test length(find_exchange_metabolites(model)["EX_fum_e"]) == 1
    
    model = convert(StandardModel, model)
    @test length(find_exchange_reactions(model)) == 21
    @test length(find_exchange_metabolites(model)["EX_fum_e"]) == 1
   
    model = convert(CoreModelCoupled, model)
    @test length(find_exchange_reactions(model)) == 21
    @test length(find_exchange_metabolites(model)["EX_fum_e"]) == 1
end
