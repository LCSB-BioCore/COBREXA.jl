@testset "Basic analysis" begin
    model = CobraTools.read_model(joinpath("data", "e_coli_core.json"))
    @test length(model.reactions) == 95 # read in correctly
    
    # FBA
    biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
    optimizer = Tulip.Optimizer # quiet by default
    sol = fba(model, biomass, optimizer)

    flux_vec = [sol[rxn.id] for rxn in model.reactions]
    sol_mapped = map_fluxes(flux_vec, model)
    @testset "FBA" begin
        @test isapprox(sol_mapped["BIOMASS_Ecoli_core_w_GAM"], 0.8739215022678488, atol=1e-6)
        @test isapprox(sol["BIOMASS_Ecoli_core_w_GAM"], 0.8739215022678488, atol=1e-6)
    end

    # exchange trackers
    atom_fluxes = atom_exchange(sol, model)
    @test isapprox(atom_fluxes["C"], -37.1902, atol=1e-3)
    
    consuming, producing = exchange_reactions(sol; verbose=false)
    @test isapprox(consuming["EX_nh4_e"], -4.76532, atol=1e-3)
    
    consuming, producing = metabolite_fluxes(sol, model)
    @test isapprox(consuming["atp_c"]["PFK"], -7.47738, atol=1e-3)
    @test isapprox(producing["atp_c"]["PYK"], 1.75818, atol=1e-3)

    # set bounds
    cbmodel, v, mb, ubs, lbs = build_cbm(model)
    glucose_index = model[findfirst(model.reactions, "EX_glc__D_e")]
    o2_index = model[findfirst(model.reactions, "EX_o2_e")]
    atpm_index = model[findfirst(model.reactions, "ATPM")]
    set_bound(glucose_index, ubs, lbs; ub=-1.0, lb=-1.0)
    @test normalized_rhs(ubs[glucose_index]) == -1.0
    @test normalized_rhs(lbs[glucose_index]) == 1.0
    set_bound(o2_index, ubs, lbs; ub=1.0, lb=1.0)
    @test normalized_rhs(ubs[o2_index]) == 1.0
    @test normalized_rhs(lbs[o2_index]) == -1.0

    # pFBA
    optimizer = OSQP.Optimizer 
    atts = Dict("eps_abs" => 5e-4,"eps_rel" => 5e-4, "max_iter" => 100_000, "verbose"=>false) 
    sol = pfba(model, biomass, optimizer; solver_attributes=atts) # just see if it works - OSQP is a terrible LP solver
    sol = pfba(model, biomass, [Tulip.Optimizer, OSQP.Optimizer]; solver_attributes=Dict("opt1" => Dict{Any, Any}(), "opt2" => atts)) # try two optimizers
    
    @testset "pFBA" begin
        @test !isempty(sol)
        @test isapprox(sol["PGM"], -14.737442319041387, atol=1e-6)        
    end


end
