@testset "Test conversion from JSONModel to StandardModel" begin

    jsonmodel = load_model(model_paths["e_coli_core.json"])
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
