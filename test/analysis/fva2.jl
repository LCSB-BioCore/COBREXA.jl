@testset "Flux variability analysis with CobraModel" begin
    model = read_model(download("http://bigg.ucsd.edu/static/models/e_coli_core.json", joinpath("data", "e_coli_core.json")))
    @test length(model.reactions) == 95 # read in correctly

    biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")

    # FVA
    optimizer = Tulip.Optimizer
    atts = Dict("IPM_IterationsLimit" => 400)
    cons = Dict("EX_glc__D_e" => (-10.0, -10.0))
    fva_max, fva_min = fva(model, biomass, optimizer, solver_attributes = atts)
    fva_max2, fva_min2 =
        fva(model, [biomass, pfl], optimizer, weights = [0.5, 0.5], constraints = cons)
    @testset "FVA" begin
        @test isapprox(fva_max["PDH"]["PDH"], 9.338922420065819, atol = 1e-6)
        @test isapprox(fva_min["PDH"]["PDH"], 9.270274952732315, atol = 1e-6)
        @test !isempty(fva_max2)
        @test !isempty(fva_min2)
    end
end
