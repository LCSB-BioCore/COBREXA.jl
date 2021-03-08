@testset "External Tests" begin
    brenda_data = parse_brenda(joinpath("data", "small_brenda.txt"))

    @test length(brenda_data) == 2
    @test brenda_data[1].ID == "1.1.1.10"
    @test brenda_data[1].TN[1].val == 0.52

    model = CobraTools.read_model(joinpath("data", "e_coli_core.json"))
    biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
    atts = Dict("eps_abs" => 5e-4,"eps_rel" => 5e-4, "max_iter" => 100_000, "verbose"=>false) 
    sol = pfba(model, biomass, [Tulip.Optimizer, OSQP.Optimizer]; solver_attributes=Dict("opt1" => Dict{Any, Any}(), "opt2" => atts)) # try two optimizers

    r_str = CobraTools.build_rxn_string(model.reactions[64])
    rs = split(replace(r_str, r"(1.0 )|( \+)|(= )|(KEGG:)"=>""), " ")
    @test all([x in ["C00354", "C00111", "C00661"] for x in rs ])
    
    gibbs_str = JSON.parsefile(joinpath("data", "gibbs.json"))
    gibbs = Dict{String, Measurement{Float64}}()
    for (k, v) in gibbs_str
        gibbs[k] = parse(Measurement{Float64}, v)
    end
    
    ΔG, mf = map_gibbs_external(sol, gibbs)
    @test isapprox(ΔG.val, -3006.0, atol=10.0)
    @test isapprox(ΔG.err, 94.0, atol=10.0)
    @test isapprox(mf, 0.034094939786929804, atol=1e-1)
        
    ΔG, mf = map_gibbs_internal(sol, gibbs)
    @test isapprox(ΔG.val, -11200.0, atol=100.0)
    @test isapprox(ΔG.err, 1300.0, atol=50.0)
    @test isapprox(mf,  0.20039554480374108, atol=1e-1)
end