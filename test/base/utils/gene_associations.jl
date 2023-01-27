
@testset "GRR parsing" begin
    @test sort(
        sort.(
            COBREXA._parse_grr("(αλφα OR βητα\x03) AND (R2-D2's_gene OR prefix:su[ff]ix)")
        ),
    ) == [
        ["R2-D2's_gene", "αλφα"],
        ["R2-D2's_gene", "βητα\x03"],
        ["prefix:su[ff]ix", "αλφα"],
        ["prefix:su[ff]ix", "βητα\x03"],
    ]
    @test_throws DomainError COBREXA._parse_grr("(")
    @test isnothing(COBREXA._parse_grr(" "))
end
