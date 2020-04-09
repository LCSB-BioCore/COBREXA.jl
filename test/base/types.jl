@testset "LinearModel type" begin
    cp = test_LP()
    @test cp isa LinearModel
end
