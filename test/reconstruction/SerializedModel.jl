
@testset "Serialized modifications" begin
    m = test_LP()

    sm = serialize_model(m, tmpfile("recon.serialized"))
    m2 = unwrap_serialized(sm)

    @test typeof(m2) == typeof(m)

    sm = serialize_model(m, tmpfile("recon.serialized"))
    m2 = remove_reaction(sm, variable_ids(m)[3])

    @test typeof(m2) == typeof(m)
    @test !(variable_ids(m)[3] in variable_ids(m2))
end
