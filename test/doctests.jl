using Test, Documenter, COBREXA

ex = quote
    using COBREXA, Tulip
    include(joinpath(dirname(pathof(COBREXA)), "..", "test", "data_static.jl"))
    model = test_LP()
end

# set module-level metadata
DocMeta.setdocmeta!(COBREXA, :DocTestSetup, ex; recursive = true)

@testset "Documentation tests" begin
    doctest(COBREXA; manual = true)
end
