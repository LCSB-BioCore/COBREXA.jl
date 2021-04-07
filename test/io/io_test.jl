@testset "IO" begin

    # E. coli models - realistic size models
    iJO1366_xml = download_data_file(
        "http://bigg.ucsd.edu/static/models/iJO1366.xml",
        joinpath("data", "iJO1366.xml"),
        "d6d9ec61ef6f155db5bb2f49549119dc13b96f6098b403ef82ea4240b27232eb",
    )
    sbmlmodel_ecoli = read_model(iJO1366_xml)
    @test_broken length(sbmlmodel_ecoli.reactions) == 2583

    iJO1366_mat = download_data_file(
        "http://bigg.ucsd.edu/static/models/iJO1366.mat",
        joinpath("data", "iJO1366.mat"),
        "b5cfe21b6369a00e45d600b783f89521f5cc953e25ee52c5f1d0a3f83743be30",
    )
    matlabmodel_ecoli = read_model(iJO1366_mat)
    @test length(matlabmodel_ecoli.reactions) == 2583

    iJO1366_json = download_data_file(
        "http://bigg.ucsd.edu/static/models/iJO1366.json",
        joinpath("data", "iJO1366.json"),
        "9376a93f62ad430719f23e612154dd94c67e0d7c9545ed9d17a4d0c347672313",
    )
    jsonmodel_ecoli = read_model(StandardModel, JSONFile, iJO1366_json)
    @test length(jsonmodel_ecoli.reactions) == 2583

    @test model_comparison_test(jsonmodel_ecoli, matlabmodel_ecoli)
    @test_broken model_comparison_test(jsonmodel_ecoli, sbmlmodel_ecoli)

    @test read_write_read_test(jsonmodel_ecoli, "json")
    @test read_write_read_test(matlabmodel_ecoli, "mat")
    @test_broken read_write_read_test(sbmlmodel_ecoli, "xml")
end
