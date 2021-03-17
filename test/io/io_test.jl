@testset "IO" begin

    # E. coli models - realistic size models
    iJO1366_xml = download(
        "http://bigg.ucsd.edu/static/models/iJO1366.xml",
        joinpath("data", "iJO1366.xml"),
    )
    sbmlmodel_ecoli = read_model(iJO1366_xml)
    @test_broken length(sbmlmodel_ecoli.reactions) == 2583

    iJO1366_mat = download(
        "http://bigg.ucsd.edu/static/models/iJO1366.mat",
        joinpath("data", "iJO1366.mat"),
    )
    matlabmodel_ecoli = read_model(iJO1366_mat)
    @test length(matlabmodel_ecoli.reactions) == 2583

    iJO1366_json = download(
        "http://bigg.ucsd.edu/static/models/iJO1366.json",
        joinpath("data", "iJO1366.json"),
    )
    jsonmodel_ecoli = read_model(iJO1366_json)
    @test length(jsonmodel_ecoli.reactions) == 2583

    @test model_comparison_test(jsonmodel_ecoli, matlabmodel_ecoli)
    @test_broken model_comparison_test(jsonmodel_ecoli, sbmlmodel_ecoli)

    @test read_write_read_test(jsonmodel_ecoli, "json")
    @test read_write_read_test(matlabmodel_ecoli, "mat")
    @test_broken read_write_read_test(sbmlmodel_ecoli, "xml")
end
