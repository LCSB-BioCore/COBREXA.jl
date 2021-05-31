
@testset "Conversion from and to MATLAB model" begin
    filename = model_paths["iJO1366.mat"]

    mm = load_mat_model(filename)
    sm = convert(StandardModel, mm)
    mm2 = convert(MATModel, sm)

    @test Set(reactions(mm)) == Set(reactions(sm))
    @test Set(reactions(mm)) == Set(reactions(mm2))
end
