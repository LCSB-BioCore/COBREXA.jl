@testset "Flux balance analysis" begin
    cp = test_simpleLP()
    (lp, x) = flux_balance_analysis(cp, GLPK.Optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(x)
    @test sol ≈ [1.0, 2.0]

    (lp, x) = flux_balance_analysis(cp, Clp.Optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(x)
    @test sol ≈ [1.0, 2.0]

    # test the maximization of the objective
    cp = test_simpleLP2()
    (lp, x) = flux_balance_analysis(cp, GLPK.Optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(x)
    @test sol ≈ [-1.0, 2.0]

    # test with a more biologically meaningfull model
    model_path = joinpath("data", "fba.mat")
    download_data_file(
        "http://bigg.ucsd.edu/static/models/iJR904.mat",
        model_path,
        "d17be86293d4caafc32b829da4e2d0d76eb45e1bb837e0138327043a83e20c6e",
    )
    cp = load_model(model_path, "iJR904")
    expected_optimum = 0.9219480950504393

    (lp, x) = flux_balance_analysis(cp, GLPK.Optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(x)
    @test objective_value(lp) ≈ expected_optimum
    @test cp.c' * sol ≈ expected_optimum

    # test the "nicer output" variants
    @test_broken false # reminder to implement these methods
    # fluxes_vec = flux_balance_analysis_vec(cp, GLPK.Optimizer)
    # @test_broken all(fluxes_vec .== sol)
    # fluxes_dict = flux_balance_analysis_dict(cp, GLPK.Optimizer)
    # rxns = reactions(cp)
    # @test all([fluxes_dict[rxns[i]] == sol[i] for i in eachindex(rxns)])
end

@testset "Flux balance analysis with CobraModel" begin
    model = read_model(
        download_data_file(
            "http://bigg.ucsd.edu/static/models/e_coli_core.json",
            joinpath("data", "e_coli_core.json"),
            "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
        ),
    )
    @test length(model.reactions) == 95 # read in correctly

    # FBA
    biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
    cons = Dict("EX_glc__D_e" => (-12.0, -12.0))
    optimizer = Tulip.Optimizer # quiet by default
    sol = fba(model, optimizer; objective_func = biomass, constraints = cons)
    pfl = findfirst(model.reactions, "PFL")
    solmulti = fba(model, optimizer; objective_func = [biomass, pfl], weights = [0.8, 0.2]) # classic flux balance analysis


    flux_vec = [sol[rxn.id] for rxn in model.reactions]
    sol_mapped = map_fluxes(flux_vec, model)
    @test isapprox(sol_mapped["BIOMASS_Ecoli_core_w_GAM"], 1.0572509997013568, atol = 1e-6)
    @test isapprox(sol["BIOMASS_Ecoli_core_w_GAM"], 1.0572509997013568, atol = 1e-6)
    @test !isempty(solmulti)

    sol = fba(model, optimizer; objective_func = biomass)

    # atom tracker
    atom_fluxes = atom_exchange(sol, model)
    @test isapprox(atom_fluxes["C"], -37.1902, atol = 1e-3)

    # exchange trackers
    consuming, producing = exchange_reactions(sol; verbose = false)
    @test isapprox(consuming["EX_nh4_e"], -4.76532, atol = 1e-3)

    # metabolite trackers
    consuming, producing = metabolite_fluxes(sol, model)
    @test isapprox(consuming["atp_c"]["PFK"], -7.47738, atol = 1e-3)
    @test isapprox(producing["atp_c"]["PYK"], 1.75818, atol = 1e-3)

    # set bounds
    cbm, v, mb, lbs, ubs = make_optimization_model(model, optimizer)
    glucose_index = model[findfirst(model.reactions, "EX_glc__D_e")]
    o2_index = model[findfirst(model.reactions, "EX_o2_e")]
    atpm_index = model[findfirst(model.reactions, "ATPM")]
    set_bound(glucose_index, lbs, ubs; ub = -1.0, lb = -1.0)
    @test normalized_rhs(ubs[glucose_index]) == -1.0
    @test normalized_rhs(lbs[glucose_index]) == 1.0
    set_bound(o2_index, lbs, ubs; ub = 1.0, lb = 1.0)
    @test normalized_rhs(ubs[o2_index]) == 1.0
    @test normalized_rhs(lbs[o2_index]) == -1.0
end
