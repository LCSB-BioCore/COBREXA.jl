using Test, Documenter, COBREXA

ex = quote
    using COBREXA, Tulip
    include(joinpath(dirname(pathof(COBREXA)), "..", "test", "data_static.jl"))
    model = test_LP()
    core_model_path = joinpath(dirname(pathof(COBREXA)), "..", "test", "data", "e_coli_core.json")
    Downloads.download(
        "http://bigg.ucsd.edu/static/models/e_coli_core.json",
        core_model_path,
    )
end

# set module-level metadata
DocMeta.setdocmeta!(COBREXA, :DocTestSetup, ex; recursive = true)

@testset "Documentation tests" begin
    doctest(COBREXA; manual = true)
end
