
@testset "Serialized models" begin
    m = load_model(CoreModel, joinpath("data", "toyModel1.mat"))

    sm = serialize_model(m, joinpath("data", "toy1.smod"))
    sm2 = serialize_model(sm, joinpath("data", "toy2.smod"))

    @test typeof(sm) == Serialized{CoreModel} # expected type
    @test typeof(sm2) == Serialized{CoreModel} # no multi-layer serialization

    precache!(sm)

    @test isequal(m, sm.m) # the data is kept okay
    @test sm2.m == nothing # nothing is cached here
    @test isequal(m, COBREXA.Serialization.deserialize(joinpath("data", "toy2.smod"))) # it was written as-is
    @test issetequal(
        reactions(convert(StandardModel, sm)),
        reactions(convert(StandardModel, sm2)),
    )
end
