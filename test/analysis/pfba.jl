@testset "Parsimonious flux balance analysis with StandardModel" begin
    model = read_model(
        download_data_file(
            "http://bigg.ucsd.edu/static/models/e_coli_core.json",
            joinpath("data", "e_coli_core.json"),
            "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
        ),
    )

    biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
    glucose = findfirst(model.reactions, "EX_glc__D_e")
    
    d = parsimonious_flux_balance_analysis_dict(model, Tulip.Optimizer; modifications=modify_constraint(glucose, -12, -12), qp_solver=modify_solver(OSQP.Optimizer), qp_solver_attributes=modify_solver_attribute("verbose", false))
    v = parsimonious_flux_balance_analysis_vec(model, Tulip.Optimizer; modifications=modify_constraint(glucose, -12, -12), qp_solver=modify_solver(OSQP.Optimizer), qp_solver_attributes=modify_solver_attribute("verbose", false))
    
    @test isapprox(d["PGM"], -14.751259785795853, atol = 1e-3)
    @test isapprox(v[8], -14.751259785795853, atol = 1e-3)    
end
