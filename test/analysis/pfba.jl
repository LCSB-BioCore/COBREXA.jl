@testset "Parsimonious flux balance analysis with CobraModel" begin

    model = read_model(
        download_data_file(
            "http://bigg.ucsd.edu/static/models/e_coli_core.json",
            joinpath("data", "e_coli_core.json"),
            "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
        ),
    )

    biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
    cons = Dict("EX_glc__D_e" => (-12.0, -12.0))

    # pFBA
    atts = Dict(
        "eps_abs" => 5e-4,
        "eps_rel" => 5e-4,
        "max_iter" => 100_000,
        "verbose" => false,
    )
    solworks = pfba(
        model,
        OSQP.Optimizer;
        objective_func = biomass,
        solver_attributes = atts,
        constraints = cons,
    ) # just see if it works - OSQP is a terrible LP solver
    sol = pfba(
        model,
        [Tulip.Optimizer, OSQP.Optimizer];
        objective_func = biomass,
        solver_attributes = Dict("opt1" => Dict{Any,Any}(), "opt2" => atts),
    ) # try two optimizers

    @test !isempty(solworks)
    @test isapprox(sol["PGM"], -14.737442319041387, atol = 1e-3)
end
