@testset "Test conversion from JSONModel to StandardModel" begin

    model_path = download_data_file(
        "http://bigg.ucsd.edu/static/models/e_coli_core.json",
        joinpath("data", "e_coli_core.json"),
        "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
    )

    jsonmodel = load_model(model_path)
    stdmodel = convert(StandardModel, jsonmodel)

    # test if same reaction ids
    @test issetequal(reactions(jsonmodel), reactions(stdmodel))
    @test issetequal(metabolites(jsonmodel), metabolites(stdmodel))
    @test issetequal(genes(jsonmodel), genes(stdmodel))
    # not the best tests since it is possible that error could cancel each other out:
    @test sum(stoichiometry(jsonmodel)) == sum(stoichiometry(stdmodel))
    jlbs, jubs = bounds(jsonmodel)
    slbs, subs = bounds(jsonmodel)
    @test sum(jlbs) == sum(slbs)
    @test sum(jubs) == sum(subs)
end
