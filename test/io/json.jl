@testset "Test conversion from JSONModel to ObjectModel" begin

    jsonmodel = load_model(model_paths["e_coli_core.json"])
    stdmodel = convert(ObjectModel, jsonmodel)

    # test if same reaction ids
    @test issetequal(variables(jsonmodel), variables(stdmodel))
    @test issetequal(metabolites(jsonmodel), metabolites(stdmodel))
    @test issetequal(genes(jsonmodel), genes(stdmodel))
    # not the best tests since it is possible that error could cancel each other out:
    @test sum(stoichiometry(jsonmodel)) == sum(stoichiometry(stdmodel))
    jlbs, jubs = bounds(jsonmodel)
    slbs, subs = bounds(jsonmodel)
    @test sum(jlbs) == sum(slbs)
    @test sum(jubs) == sum(subs)
end

@testset "Save JSON model" begin
    model = load_model(MatrixModel, model_paths["e_coli_core.json"])
    testpath = tmpfile("modeltest.json")
    save_model(model, testpath)
    wrote = convert(MatrixModel, load_json_model(testpath))
    @test isequal(model, wrote)
end
