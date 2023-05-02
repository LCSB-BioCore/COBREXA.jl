
struct FakeModel <: AbstractMetabolicModel
    dummy::Int
end

@testset "Base abstract model methods require proper minimal implementation" begin
    @test_throws MethodError variable_ids(123)
    x = FakeModel(123)
    for m in [variables, metabolites, stoichiometry, bounds, objective]
        @test_throws MethodError m(x)
    end
end


@testset "ID shortcuts are identictical with the ID-generating functions" begin
    @test variables === variable_ids
    @test reactions === reaction_ids
    # TODO don't forget about metabolites later
end
