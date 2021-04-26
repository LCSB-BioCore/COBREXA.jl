@testset "Opening models from BIGG" begin

    # E. coli models - realistic size models
    iJO1366_xml = download_data_file(
        "http://bigg.ucsd.edu/static/models/iJO1366.xml",
        joinpath("data", "iJO1366.xml"),
        "d6d9ec61ef6f155db5bb2f49549119dc13b96f6098b403ef82ea4240b27232eb",
    )
    sbmlmodel = load_model(iJO1366_xml)
    @test sbmlmodel isa SBMLModel
    @test n_reactions(sbmlmodel) == 2583

    iJO1366_mat = download_data_file(
        "http://bigg.ucsd.edu/static/models/iJO1366.mat",
        joinpath("data", "iJO1366.mat"),
        "b5cfe21b6369a00e45d600b783f89521f5cc953e25ee52c5f1d0a3f83743be30",
    )
    matlabmodel = load_model(iJO1366_mat)
    @test matlabmodel isa MATModel
    @test n_reactions(matlabmodel) == 2583

    iJO1366_json = download_data_file(
        "http://bigg.ucsd.edu/static/models/iJO1366.json",
        joinpath("data", "iJO1366.json"),
        "9376a93f62ad430719f23e612154dd94c67e0d7c9545ed9d17a4d0c347672313",
    )
    jsonmodel = load_model(iJO1366_json)
    @test jsonmodel isa JSONModel
    @test n_reactions(jsonmodel) == 2583

    @test Set(lowercase.(reactions(sbmlmodel))) ==
          Set("r_" .* lowercase.(reactions(matlabmodel)))
    @test Set(lowercase.(reactions(sbmlmodel))) ==
          Set("r_" .* lowercase.(reactions(jsonmodel)))
end
