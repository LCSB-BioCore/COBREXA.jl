@testset "CoreModel type" begin
    cp = test_LP()
    @test cp isa CoreModel
end
