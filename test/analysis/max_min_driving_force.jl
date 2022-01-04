@testset "Max-min driving force analysis" begin

    model = load_model(model_paths["e_coli_core.json"])

    flux_solution = flux_balance_analysis_dict(
        model,
        GLPK.Optimizer;
        modifications = [add_loopless_constraints()],
    )

    sol = max_min_driving_force(
        model,
        reaction_standard_gibbs_free_energies,
        Tulip.Optimizer;
        flux_solution = flux_solution,
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
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    @test isapprox(sol.mmdf, 1.7661155558545698, atol = TEST_TOLERANCE)

    sols = max_min_driving_force_variability(
        model,
        reaction_standard_gibbs_free_energies,
        Tulip.Optimizer;
        bounds = gamma_bounds(0.9),
        flux_solution = flux_solution,
        proton_ids = ["h_c", "h_e"],
        water_ids = ["h2o_c", "h2o_e"],
        concentration_ratios = Dict{Tuple{String,String},Float64}(
            ("atp_c", "adp_c") => 10.0,
            ("nadh_c", "nad_c") => 0.13,
            ("nadph_c", "nadp_c") => 1.3,
        ),
        constant_concentrations = Dict{String,Float64}(
        # "pi_c" => 10e-3
        ),
        concentration_lb = 1e-6,
        concentration_ub = 100e-3,
        ignore_reaction_ids = ["H2Ot"],
        modifications = [change_optimizer_attribute("IPM_IterationsLimit", 1000)],
    )

    pyk_idx = first(indexin(["PYK"], reactions(model)))
    @test isapprox(
        sols[pyk_idx, 1].dg_reactions["PYK"],
        -1.5895040002691128;
        atol = TEST_TOLERANCE,
    )
end
