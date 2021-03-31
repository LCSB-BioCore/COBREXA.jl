@testset "Flux balance analysis" begin
    cp = test_simpleLP()
    (lp, x) = fluxBalanceAnalysis(cp, GLPK.Optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(x)
    @test sol ≈ [1.0, 2.0]

    (lp, x) = fluxBalanceAnalysis(cp, Clp.Optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(x)
    @test sol ≈ [1.0, 2.0]

    # test the maximization of the objective
    cp = test_simpleLP2()
    (lp, x) = fluxBalanceAnalysis(cp, GLPK.Optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(x)
    @test sol ≈ [-1.0, 2.0]

    # test with a more biologically meaningfull model
    modelPath = joinpath("data", "fba.mat")
    downloadDataFile(
        "http://bigg.ucsd.edu/static/models/iJR904.mat",
        modelPath,
        "d17be86293d4caafc32b829da4e2d0d76eb45e1bb837e0138327043a83e20c6e",
    )
    cp = loadModel(modelPath, "iJR904")
    expectedOptimum = 0.9219480950504393

    (lp, x) = fluxBalanceAnalysis(cp, GLPK.Optimizer)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(x)
    @test objective_value(lp) ≈ expectedOptimum
    @test cp.c' * sol ≈ expectedOptimum

    # test the "nicer output" variants
    @test_broken false # reminder to implement these methods
    # fluxesVec = fluxBalanceAnalysisVec(cp, GLPK.Optimizer)
    # @test_broken all(fluxesVec .== sol)
    # fluxesDict = fluxBalanceAnalysisDict(cp, GLPK.Optimizer)
    # rxns = reactions(cp)
    # @test all([fluxesDict[rxns[i]] == sol[i] for i in eachindex(rxns)])
end

@testset "Flux balance analysis with CobraModel" begin
    model = read_model(
        Downloads.download(
            "http://bigg.ucsd.edu/static/models/e_coli_core.json",
            joinpath("data", "e_coli_core.json"),
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
    cbm, v, mb, lbs, ubs = makeOptimizationModel(model, optimizer)
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
