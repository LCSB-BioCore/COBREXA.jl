
@testset "Serialized models" begin
    m = test_simpleLP()

    sm = serialize_model(m, tmpfile("toy1.serialized"))
    sm2 = serialize_model(sm, tmpfile("toy2.serialized"))

    @test typeof(sm) == Serialized{CoreModel} # expected type
    @test typeof(sm2) == Serialized{CoreModel} # no multi-layer serialization

    precache!(sm)

    @test isequal(m, sm.m) # the data is kept okay
    @test sm2.m == nothing # nothing is cached here
    @test isequal(m, COBREXA.Serialization.deserialize(tmpfile("toy2.serialized"))) # it was written as-is
    @test issetequal(
        reactions(convert(StandardModel, sm)),
        reactions(convert(StandardModel, sm2)),
    )
    sm.m = nothing
    @test issetequal(
        metabolites(convert(CoreModelCoupled, sm)),
        metabolites(convert(CoreModelCoupled, sm2)),
    )
end
