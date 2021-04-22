@testset "Flux balance analysis with CoreModel" begin
    cp = test_simpleLP()
    lp = flux_balance_analysis(cp, GLPK.Optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(lp[:x])
    @test sol ≈ [1.0, 2.0]

    lp = flux_balance_analysis(cp, Clp.Optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(lp[:x])
    @test sol ≈ [1.0, 2.0]

    # test the maximization of the objective
    cp = test_simpleLP2()
    lp = flux_balance_analysis(cp, GLPK.Optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(lp[:x])
    @test sol ≈ [-1.0, 2.0]

    # test with a more biologically meaningfull model
    model_path = joinpath("data", "fba.mat")
    download_data_file(
        "http://bigg.ucsd.edu/static/models/iJR904.mat",
        model_path,
        "d17be86293d4caafc32b829da4e2d0d76eb45e1bb837e0138327043a83e20c6e",
    )
    cp = read_model(model_path, CoreModel)
    expected_optimum = 0.9219480950504393

    lp = flux_balance_analysis(cp, GLPK.Optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(lp[:x])
    @test objective_value(lp) ≈ expected_optimum
    @test cp.c' * sol ≈ expected_optimum

    # test the "nicer output" variants
    fluxes_vec = flux_balance_analysis_vec(cp, GLPK.Optimizer)
    @test all(fluxes_vec .== sol)
    fluxes_dict = flux_balance_analysis_dict(cp, GLPK.Optimizer)
    rxns = reactions(cp)
    @test all([fluxes_dict[rxns[i]] == sol[i] for i in eachindex(rxns)])
end

@testset "Flux balance analysis with StandardModel" begin
    model_path = download_data_file(
        "http://bigg.ucsd.edu/static/models/e_coli_core.json",
        joinpath("data", "e_coli_core.json"),
        "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
    )

    model = read_model(model_path, StandardModel)

    biomass = model.reactions["BIOMASS_Ecoli_core_w_GAM"]
    glucose = model.reactions["EX_glc__D_e"]
    sol = flux_balance_analysis_dict(
        model,
        Tulip.Optimizer;
        modifications = [
            change_objective(biomass),
            change_constraint(glucose, -12, -12),
            change_sense(MOI.MAX_SENSE),
            change_solver_attribute("IPM_IterationsLimit", 110),
        ],
    )
    @test isapprox(
        sol["BIOMASS_Ecoli_core_w_GAM"],
        1.0572509997013568,
        atol = TEST_TOLERANCE,
    )

    pfl = model.reactions["PFL"]
    pfl_frac = 0.8
    biomass_frac = 0.2
    sol_multi = flux_balance_analysis_dict(
        model,
        Tulip.Optimizer;
        modifications = change_objective(
            [biomass, pfl];
            weights = [biomass_frac, pfl_frac],
        ),
    )
    @test isapprox(
        biomass_frac * sol_multi["BIOMASS_Ecoli_core_w_GAM"] + pfl_frac * sol_multi["PFL"],
        31.999999998962604,
        atol = TEST_TOLERANCE,
    )
end
