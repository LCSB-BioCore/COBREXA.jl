@testset "Brenda Tests" begin
    brenda_data = parse_brenda(joinpath("data", "small_brenda.txt"))

    @test length(brenda_data) == 2
    @test brenda_data[1].ID == "1.1.1.10"
    @test brenda_data[1].TN[1].val == 0.52
end