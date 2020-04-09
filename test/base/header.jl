
@testset "Header" begin
    loadRes = COBREXA.loadSource(["base", "io", "reconstruction", "analysis"], joinpath(dirname(pathof(COBREXA))), true)
    @test any(val->!val, loadRes) == false
end