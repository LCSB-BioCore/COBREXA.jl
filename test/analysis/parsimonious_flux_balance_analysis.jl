@testset "Parsimonious flux balance analysis with StandardModel" begin
    model_path = download_data_file(
        "http://bigg.ucsd.edu/static/models/e_coli_core.json",
        joinpath("data", "e_coli_core.json"),
        "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
    )

    model = read_model(model_path, StandardModel)

    biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
    glucose = findfirst(model.reactions, "EX_glc__D_e")

    d = parsimonious_flux_balance_analysis_dict(
        model,
        Tulip.Optimizer;
        modifications = [
            change_constraint(glucose, -12, -12),
            change_solver_attribute("IPM_IterationsLimit", 500),
        ],
        qp_solver = change_solver(OSQP.Optimizer),
        qp_solver_attributes = change_solver_attribute("verbose", false),
    )
    v = parsimonious_flux_balance_analysis_vec(
        model,
        Tulip.Optimizer;
        modifications = [
            change_constraint(glucose, -12, -12),
            change_solver_attribute("IPM_IterationsLimit", 500),
        ],
        qp_solver = change_solver(OSQP.Optimizer),
        qp_solver_attributes = change_solver_attribute("verbose", false),
    )

    @test isapprox(d["PGM"], -17.568590034769613, atol = 0.5)
    @test isapprox(v[8], -17.568590034769613, atol = 0.5) # OSQP sucks
end
