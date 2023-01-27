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

    # specifically test parsing of gene-reaction associations in Recon
    reconmodel = load_model(StandardModel, model_paths["Recon3D.json"])
    @test n_reactions(reconmodel) == 10600
    recon_grrs = [r.grr for (i, r) in reconmodel.reactions if !isnothing(r.grr)]
    @test length(recon_grrs) == 5938
    @test sum(length.(recon_grrs)) == 31504
end
