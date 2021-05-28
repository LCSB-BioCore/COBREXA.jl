agora_model = joinpath("data", "agora-model.mat")

download_data_file(
    "https://www.vmh.life/files/reconstructions/AGORA/1.03/reconstructions/mat/Mycoplasma_hominis_ATCC_23114.mat",
    agora_model,
    "03362073aa917f0691a0c896948f6e8eebe47f02dcbe0c3f00275fa87396e220",
)

@testset "Import MAT model" begin
    cp = load_model(CoreModel, agora_model)
    @test cp isa CoreModel
    @test size(cp.S) == (475, 496)
end

@testset "Save MAT model" begin
    testpath = joinpath("data", "agora-clone.mat")
    loaded = load_model(CoreModel, agora_model)
    save_model(loaded, testpath)
    wrote = load_model(CoreModel, testpath)
    @test wrote isa CoreModel
    @test isequal(wrote, loaded)
end
