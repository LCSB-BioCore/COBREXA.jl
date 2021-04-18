
struct FakeModel <: MetabolicModel
    dummy::Int
end

@testset "Base abstract model methods require proper implementation" begin
    @test_throws MethodError reactions(123)
    x = FakeModel(123)
    for m in [reactions, metabolites, stoichiometry, bounds, balance, objective]
        @test_throws MethodError m(x)
    end
end
