using Test, Documenter, COBREXA

ex = quote
    using COBREXA, Tulip
    include(joinpath(dirname(pathof(COBREXA)), "..", "test", "data_static.jl"))
    model = test_LP()
    include(joinpath(dirname(pathof(COBREXA)), "..", "test", "test_functions.jl")) # expose functions for downloading
    include(joinpath(dirname(pathof(COBREXA)), "..", "test", "doctest_models.jl")) 
end

# set module-level metadata
DocMeta.setdocmeta!(COBREXA, :DocTestSetup, ex; recursive = true)

@testset "Documentation tests" begin
    doctest(COBREXA; manual = true)
end
