
@testset "Serialized models" begin
    m = test_simpleLP()

    sm = serialize_model(m, tmpfile("toy1.serialized"))
    sm2 = serialize_model(sm, tmpfile("toy2.serialized"))

    @test typeof(sm) == Serialized{MatrixModel} # expected type
    @test typeof(sm2) == Serialized{MatrixModel} # no multi-layer serialization

    precache!(sm)

    @test isequal(m, sm.m) # the data is kept okay
    @test sm2.m == nothing # nothing is cached here
    @test isequal(m, deserialize(tmpfile("toy2.serialized"))) # it was written as-is
    @test issetequal(
        variable_ids(convert(ObjectModel, sm)),
        variable_ids(convert(ObjectModel, sm2)),
    )
    sm.m = nothing
    @test issetequal(
        metabolites(convert(MatrixModelWithCoupling, sm)),
        metabolites(convert(MatrixModelWithCoupling, sm2)),
    )
end
