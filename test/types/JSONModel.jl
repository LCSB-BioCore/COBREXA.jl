@testset "Conversion from and to SBML model" begin
    json_model = model_paths["iJO1366.json"]

    jm = load_json_model(json_model)
    sm = convert(ObjectModel, jm)
    jm2 = convert(JSONModel, sm)

    @test Set(variables(jm)) == Set(variables(sm))
    @test Set(variables(jm)) == Set(variables(jm2))
end

@testset "JSONModel generic interface" begin
    model = load_model(model_paths["e_coli_core.json"])

    @test reaction_stoichiometry(model, "EX_ac_e") == Dict("ac_e" => -1)
end