@testset "Flux balance analysis with CoreModel" begin
    cp = test_simpleLP()
    lp = flux_balance_analysis(cp, Tulip.Optimizer)
    @test termination_status(lp) == MOI.OPTIMAL
    sol = JuMP.value.(lp[:x])
    @test sol ≈ [1.0, 2.0]

    # test the maximization of the objective
    cp = test_simpleLP2()
    lp = flux_balance_analysis(cp, Tulip.Optimizer)
    @test termination_status(lp) == MOI.OPTIMAL
    sol = JuMP.value.(lp[:x])
    @test sol ≈ [-1.0, 2.0]

    # test with a more biologically meaningfull model
    cp = load_model(CoreModel, model_paths["iJR904.mat"])
    expected_optimum = 0.9219480950504393

    lp = flux_balance_analysis(cp, Tulip.Optimizer)
    @test termination_status(lp) == MOI.OPTIMAL
    sol = JuMP.value.(lp[:x])
    @test isapprox(objective_value(lp), expected_optimum, atol = TEST_TOLERANCE)
    @test isapprox(cp.c' * sol, expected_optimum, atol = TEST_TOLERANCE)

    # test the "nicer output" variants
    fluxes_vec = flux_balance_analysis_vec(cp, Tulip.Optimizer)
    @test all(fluxes_vec .== sol)
    fluxes_dict = flux_balance_analysis_dict(cp, Tulip.Optimizer)
    rxns = reactions(cp)
    @test all([fluxes_dict[rxns[i]] == sol[i] for i in eachindex(rxns)])
end

@testset "Flux balance analysis with StandardModel" begin

    model = load_model(StandardModel, model_paths["e_coli_core.json"])

    sol = flux_balance_analysis_dict(
        model,
        Tulip.Optimizer;
        modifications = [
            change_objective("BIOMASS_Ecoli_core_w_GAM"),
            change_constraint("EX_glc__D_e", -12, -12),
            change_sense(MAX_SENSE),
            change_optimizer_attribute("IPM_IterationsLimit", 110),
        ],
    )
    @test isapprox(
        sol["BIOMASS_Ecoli_core_w_GAM"],
        1.0572509997013568,
        atol = TEST_TOLERANCE,
    )

    pfl_frac = 0.8
    biomass_frac = 0.2
    sol_multi = flux_balance_analysis_dict(
        model,
        Tulip.Optimizer;
        modifications = [
            change_objective(
                ["BIOMASS_Ecoli_core_w_GAM", "PFL"];
                weights = [biomass_frac, pfl_frac],
            ),
        ],
    )
    @test isapprox(
        biomass_frac * sol_multi["BIOMASS_Ecoli_core_w_GAM"] + pfl_frac * sol_multi["PFL"],
        31.999999998962604,
        atol = TEST_TOLERANCE,
    )

    @test_throws DomainError flux_balance_analysis_dict(
        model,
        Tulip.Optimizer;
        modifications = [change_constraint("gbbrsh", -12, -12)],
    )
    @test_throws DomainError flux_balance_analysis_dict(
        model,
        Tulip.Optimizer;
        modifications = [change_objective("gbbrsh")],
    )
    @test_throws DomainError flux_balance_analysis_dict(
        model,
        Tulip.Optimizer;
        modifications = [change_objective(["BIOMASS_Ecoli_core_w_GAM"; "gbbrsh"])],
    )
end

@testset "Flux balance analysis with CoreModelCoupled" begin

    model = load_model(CoreModel, model_paths["e_coli_core.json"])

    # assume coupling constraints of the form:
    # -γ ≤ vᵢ/μ  ≤ γ
    # I.e., enforces that the ratio between any reaction flux
    # and the growth rate is bounded by γ.
    γ = 40

    # construct coupling bounds
    nr = n_reactions(model)
    biomass_index = first(indexin(["BIOMASS_Ecoli_core_w_GAM"], reactions(model)))

    Cf = sparse(1.0I, nr, nr)
    Cf[:, biomass_index] .= -γ

    Cb = sparse(1.0I, nr, nr)
    Cb[:, biomass_index] .= γ

    C = [Cf; Cb]

    clb = spzeros(2 * nr)
    clb[1:nr] .= -1000.0
    cub = spzeros(2 * nr)
    cub[nr+1:end] .= 1000

    cmodel = CoreModelCoupled(model, C, clb, cub) # construct

    dc = flux_balance_analysis_dict(cmodel, Tulip.Optimizer)
    @test isapprox(dc["BIOMASS_Ecoli_core_w_GAM"], 0.665585699298256, atol = TEST_TOLERANCE)
end
