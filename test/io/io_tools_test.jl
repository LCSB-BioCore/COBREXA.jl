@testset "IO" begin

    # E. coli models - realistic size model
    iJO1366_xml = joinpath("data", "iJO1366.xml")
    sbmlmodel_ecoli = CobraTools.read_model(iJO1366_xml)
    @test_broken length(sbmlmodel_ecoli.reactions) == 2583

    iJO1366_mat = joinpath("data", "iJO1366.mat")
    matlabmodel_ecoli = CobraTools.read_model(iJO1366_mat)
    @test length(matlabmodel_ecoli.reactions) == 2583

    iJO1366_json = joinpath("data", "iJO1366.json")
    jsonmodel_ecoli = CobraTools.read_model(iJO1366_json)
    @test length(jsonmodel_ecoli.reactions) == 2583

    @test model_comparison_test(jsonmodel_ecoli, matlabmodel_ecoli)
    @test_broken model_comparison_test(jsonmodel_ecoli, sbmlmodel_ecoli)

    @test read_write_read_test(jsonmodel_ecoli, "json")
    @test read_write_read_test(matlabmodel_ecoli, "mat")
    @test_broken read_write_read_test(sbmlmodel_ecoli, "xml")
end
