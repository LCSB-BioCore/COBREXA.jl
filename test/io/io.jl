@testset "Opening models from BIGG" begin

    sbmlmodel = load_model(model_paths["iJO1366.xml"])
    @test sbmlmodel isa SBMLModel
    @test n_reactions(sbmlmodel) == 2583

    matlabmodel = load_model(model_paths["iJO1366.mat"])
    @test matlabmodel isa MATModel
    @test n_reactions(matlabmodel) == 2583

    jsonmodel = load_model(model_paths["iJO1366.json"])
    @test jsonmodel isa JSONModel
    @test n_reactions(jsonmodel) == 2583

    @test Set(lowercase.(reactions(sbmlmodel))) ==
          Set("r_" .* lowercase.(reactions(matlabmodel)))
    @test Set(lowercase.(reactions(sbmlmodel))) ==
          Set("r_" .* lowercase.(reactions(jsonmodel)))
end
