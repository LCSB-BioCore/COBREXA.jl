
@testset "Conversion from and to MATLAB model" begin
    filename = download_data_file(
        "http://bigg.ucsd.edu/static/models/iJO1366.mat",
        joinpath("data", "iJO1366.mat"),
        "b5cfe21b6369a00e45d600b783f89521f5cc953e25ee52c5f1d0a3f83743be30",
    )

    mm = load_mat_model(filename)
    cm = convert(CoreModel, mm)
    mm2 = convert(MATModel, cm)

    @test Set(reactions(mm)) == Set(reactions(mm2))
end
