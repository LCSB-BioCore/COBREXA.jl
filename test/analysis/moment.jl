@testset "Moment algorithm" begin
    model = load_model(StandardModel, model_paths["e_coli_core.json"])

    ksas = Dict(rid => 1000.0 for rid in reactions(model))
    protein_mass_fraction = 0.56

    sol = flux_balance_analysis_dict(
        model,
        Tulip.Optimizer;
        modifications = [
            add_moment_constraints(ksas, protein_mass_fraction;),
            change_constraint("EX_glc__D_e", lb = -1000),
        ],
    )

    @test isapprox(
        sol["BIOMASS_Ecoli_core_w_GAM"],
        0.6623459899423948,
        atol = TEST_TOLERANCE,
    )
end
