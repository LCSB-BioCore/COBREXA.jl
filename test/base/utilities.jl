@testset "CoreModel simple functions" begin
    cp = test_LP()
    @test n_reactions(cp) == 3
    @test n_metabolites(cp) == 4
    @test n_coupling_constraints(cp) == 0

    cp2 = test_LP()
    @test isequal(cp, cp2)
    cp2.S[1] = 1
    @test !isequal(cp, cp2)
    @test isequal(cp, copy(cp))

    cp = test_coupledLP()
    @test n_coupling_constraints(cp) == 2000
    @test isequal(cp, copy(cp))
end

@testset "Analysis utilities" begin
    model_path = download_data_file(
        "http://bigg.ucsd.edu/static/models/e_coli_core.json",
        joinpath("data", "e_coli_core.json"),
        "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
    )

    model = read_model(model_path, StandardModel)

    # FBA
    biomass = model.reactions["BIOMASS_Ecoli_core_w_GAM"]

    glc = model.reactions["EX_glc__D_e"]
    optimizer = Tulip.Optimizer # quiet by default
    sol = flux_balance_analysis_dict(
        model,
        optimizer;
        modifications = [change_objective(biomass)],
    )

    # atom tracker
    atom_fluxes = atom_exchange(sol, model)
    @test isapprox(atom_fluxes["C"], -37.1902, atol = 1e-3)

    # exchange trackers
    consuming, producing = exchange_reactions(sol, model; verbose = false)
    @test isapprox(consuming["EX_nh4_e"], -4.76532, atol = 1e-3)

    # metabolite trackers
    consuming, producing = metabolite_fluxes(sol, model)
    @test isapprox(consuming["atp_c"]["PFK"], -7.47738, atol = 1e-3)
    @test isapprox(producing["atp_c"]["PYK"], 1.75818, atol = 1e-3)

    # set bounds
    cbm = make_optimization_model(model, optimizer)
    ubs = cbm[:ubs]
    lbs = cbm[:lbs]
    glucose_index = index_of("EX_glc__D_e", reactions(model))
    o2_index = index_of("EX_o2_e", reactions(model))
    atpm_index = index_of("ATPM", reactions(model))
    set_bound(glucose_index, cbm; ub = -1.0, lb = -1.0)
    @test normalized_rhs(ubs[glucose_index]) == -1.0
    @test normalized_rhs(lbs[glucose_index]) == 1.0
    set_bound(o2_index, cbm; ub = 1.0, lb = 1.0)
    @test normalized_rhs(ubs[o2_index]) == 1.0
    @test normalized_rhs(lbs[o2_index]) == -1.0
end
