@testset "Reaction utilities" begin
    # test if reaction equation can be built back into a sensible reaction string
    req = Dict("coa_c" => -1, "for_c" => 1, "accoa_c" => 1, "pyr_c" => -1)
    rstr_out = stoichiometry_string(req)
    @test occursin("coa_c", split(rstr_out, " = ")[1])
    @test occursin("for", split(rstr_out, " = ")[2])
end
