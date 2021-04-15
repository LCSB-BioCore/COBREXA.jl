using Test, Documenter, COBREXA
@testset "Documentation tests" begin
    doctest(COBREXA; manual = false)
end
