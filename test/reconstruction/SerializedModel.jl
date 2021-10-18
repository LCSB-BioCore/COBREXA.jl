
@testset "Serialized modifications" begin
    m = test_LP()

    sm = serialize_model(m, tmpfile("recon.serialized"))
    m2 = unwrap_serialized(sm)

    @test typeof(m2) == typeof(m)

    sm = serialize_model(m, tmpfile("recon.serialized"))
    m2 = remove_reaction(sm, reactions(m)[3])

    @test typeof(m2) == typeof(m)
    @test !(reactions(m)[3] in reactions(m2))
end
