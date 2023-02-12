@testset "Max-min driving force analysis" begin

    model = load_model(model_paths["e_coli_core.json"])

    flux_solution = flux_balance_analysis_dict(
        model,
        GLPK.Optimizer;
        modifications = [add_loopless_constraints()],
    )

    mmdfm = make_max_min_driving_force_model(
        model;
        reaction_standard_gibbs_free_energies,
        flux_solution,
        proton_ids = ["h_c", "h_e"],
        water_ids = ["h2o_c", "h2o_e"],
        concentration_ratios = Dict(
            ("atp_c", "adp_c") => 10.0,
            ("nadh_c", "nad_c") => 0.13,
            ("nadph_c", "nadp_c") => 1.3,
        ),
        concentration_lb = 1e-6,
        concentration_ub = 100e-3,
        ignore_reaction_ids = ["H2Ot"],
    )

    opt_model = flux_balance_analysis(
        mmdfm,
        Tulip.Optimizer;
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    # get mmdf
    @test isapprox(
        solved_objective_value(opt_model),
        1.7661155558545698,
        atol = TEST_TOLERANCE,
    )

    # values_dict(:reaction, mmdfm, opt_model) # TODO throw missing semantics error
    @test length(values_dict(:metabolite_log_concentration, mmdfm, opt_model)) == 72
    @test length(values_dict(:gibbs_free_energy_reaction, mmdfm, opt_model)) == 95

    sols = variability_analysis(
        Val(:gibbs_free_energy_reaction),
        mmdfm,
        Tulip.Optimizer;
        bounds = gamma_bounds(0.9),
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    pyk_idx = first(indexin(["ΔG PYK"], gibbs_free_energy_reactions(mmdfm)))
    @test isapprox(sols[pyk_idx, 2], -1.5895040002691128; atol = TEST_TOLERANCE)
end
