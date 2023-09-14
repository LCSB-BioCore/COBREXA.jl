@testset "Opening models from BIGG" begin

    sbmlmodel = load_model(model_paths["iJO1366.xml"])
    @test sbmlmodel isa SBMLModel
    @test variable_count(sbmlmodel) == 2583

    matlabmodel = load_model(model_paths["iJO1366.mat"])
    @test matlabmodel isa MATModel
    @test variable_count(matlabmodel) == 2583

    jsonmodel = load_model(model_paths["iJO1366.json"])
    @test jsonmodel isa JSONModel
    @test variable_count(jsonmodel) == 2583

    @test Set(lowercase.(variable_ids(sbmlmodel))) ==
          Set("r_" .* lowercase.(variable_ids(matlabmodel)))
    @test Set(lowercase.(variable_ids(sbmlmodel))) ==
          Set("r_" .* lowercase.(variable_ids(jsonmodel)))

    # specifically test parsing of gene-reaction associations in Recon
    reconmodel = load_model(ObjectModel, model_paths["Recon3D.json"])
    @test variable_count(reconmodel) == 10600
    recon_grrs = [
        r.gene_associations for
        (i, r) in reconmodel.reactions if !isnothing(r.gene_associations)
    ]
    @test length(recon_grrs) == 5938
    @test sum(length.(recon_grrs)) == 31504
end
