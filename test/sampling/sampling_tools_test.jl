@testset "Sampling Tests" begin
    # these tests are not very good - sampling needs work
    model = read_model(download("http://bigg.ucsd.edu/static/models/e_coli_core.json", joinpath("data", "e_coli_core.json")))
    @test length(model.reactions) == 95 # read in correctly


    optimizer = Tulip.Optimizer
    biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
    cons = Dict("EX_glc__D_e" => (-12.0, -12.0))
    atts = Dict("IPM_IterationsLimit" => 110)
    sol = fba(model, optimizer; objective_func=biomass, constraints = cons)
    cons["BIOMASS_Ecoli_core_w_GAM"] =
        (sol["BIOMASS_Ecoli_core_w_GAM"] * 0.99, sol["BIOMASS_Ecoli_core_w_GAM"])
    
    samples = hit_and_run(
        100_000,
        model,
        optimizer;
        keepevery = 10,
        samplesize = 5000,
        constraints = cons,
        solver_attributes = atts,
    )
    
    # @test isapprox(mean(samples[64, :]), 8.9, atol = 0.5) # only tests if the sampler approximately converged
    
    # samples = achr(
    #     100_000,
    #     model,
    #     optimizer;
    #     keepevery = 10,
    #     samplesize = 5000,
    #     constraints = cons,
    # )
    # @test isapprox(mean(samples[64, :]), 8.9, atol = 0.5) # only tests if the sampler approximately converged

end
