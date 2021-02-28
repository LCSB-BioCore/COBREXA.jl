@testset "Basic analysis" begin
    model = CobraTools.read_model(joinpath("data", "e_coli_core.json"))
    @test length(model.reactions) == 95 # read in correctly
    
    biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
    optimizer = Tulip.Optimizer # quiet by default
    sol = fba(model, biomass, optimizer)
    @test isapprox(sol["BIOMASS_Ecoli_core_w_GAM"], 0.8739215022678488, atol=1e-6)

    optimizer = OSQP.Optimizer 
    atts = Dict("eps_abs" => 5e-4,"eps_rel" => 5e-4, "max_iter" => 100_000, "verbose"=>false) 
    sol = pfba(model, biomass, optimizer; solver_attributes=atts) # just see if it works - OSQP is a terrible LP solver
    @test !isempty(sol)

    sol = pfba(model, biomass, [Tulip.Optimizer, OSQP.Optimizer]; solver_attributes=Dict("opt1" => Dict{Any, Any}(), "opt2" => atts)) # try two optimizers
    @test isapprox(sol["PGM"], -14.737442319041387, atol=1e-6)
end
