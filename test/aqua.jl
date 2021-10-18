
@testset "Automated QUality Assurance" begin
    # We can't do Aqua.test_all here (yet) because the ambiguity tests fail in
    # deps. Instead let's pick the other sensible tests.
    Aqua.test_deps_compat(COBREXA)
    Aqua.test_project_extras(COBREXA)
    Aqua.test_project_toml_formatting(COBREXA)
    #Aqua.test_stale_deps(COBREXA) # currently seems broken
    Aqua.test_unbound_args(COBREXA)
    Aqua.test_undefined_exports(COBREXA)
end
