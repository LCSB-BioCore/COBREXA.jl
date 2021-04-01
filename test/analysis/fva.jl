@testset "Flux variability analysis" begin
    cp = test_simpleLP()
    optimizer = GLPK.Optimizer
    fluxes = flux_variability_analysis(cp, optimizer)

    @test size(fluxes) == (2, 2)
    @test fluxes ≈ [
        1.0 1.0
        2.0 2.0
    ]

    fluxes = flux_variability_analysis(cp, [2], optimizer)

    @test size(fluxes) == (1, 2)
    @test fluxes == Array{Float64,2}([2 2])

    # a special testcase for slightly sub-optimal FVA (gamma<1)
    cp = LinearModel(
        [-1.0 -1.0 -1.0],
        [0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, -1.0],
        1.0 * ones(3),
        ["r$x" for x = 1:3],
        ["m1"],
    )
    fluxes = flux_variability_analysis(cp, optimizer)
    @test fluxes ≈ [
        1.0 1.0
        0.0 0.0
        -1.0 -1.0
    ]
    fluxes = flux_variability_analysis(cp, optimizer; gamma = 0.5)
    @test fluxes ≈ [
        0.5 1.0
        0.0 0.5
        -1.0 -0.5
    ]
    fluxes = flux_variability_analysis(cp, optimizer; gamma = 0.0)
    @test fluxes ≈ [
        0.0 1.0
        0.0 1.0
        -1.0 0.0
    ]

    @test isempty(flux_variability_analysis(cp, Vector{Int}(), GLPK.Optimizer))
    @test_throws DomainError flux_variability_analysis(cp, [-1], GLPK.Optimizer)
    @test_throws DomainError flux_variability_analysis(cp, [99999999], GLPK.Optimizer)
end

@testset "Parallel FVA" begin
    cp = test_simpleLP()
    pids = addprocs(2, topology = :master_worker)
    @everywhere using COBREXA, GLPK
    fluxes = flux_variability_analysis(cp, [1, 2], GLPK.Optimizer, pids)
    @test fluxes ≈ [
        1.0 1.0
        2.0 2.0
    ]
    rmprocs(pids)
end

@testset "Flux variability analysis with CobraModel" begin
    model = read_model(
        download_data_file(
            "http://bigg.ucsd.edu/static/models/e_coli_core.json",
            joinpath("data", "e_coli_core.json"),
            "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
        ),
    )
    @test length(model.reactions) == 95 # read in correctly

    biomass = findfirst(model.reactions, "BIOMASS_Ecoli_core_w_GAM")
    pfl = findfirst(model.reactions, "PFL")

    # FVA
    optimizer = Tulip.Optimizer
    atts = Dict("IPM_IterationsLimit" => 500)
    cons = Dict("EX_glc__D_e" => (-10.0, -10.0))
    fva_max, fva_min =
        fva(model, optimizer; objective_func = biomass, solver_attributes = atts)
    fva_max2, fva_min2 = fva(
        model,
        optimizer;
        objective_func = [biomass, pfl],
        weights = [0.5, 0.5],
        constraints = cons,
    )
    @testset "FVA" begin
        @test isapprox(fva_max["PDH"]["PDH"], 9.338922420065819, atol = 1e-6)
        @test isapprox(fva_min["PDH"]["PDH"], 9.270274952732315, atol = 1e-6)
        @test !isempty(fva_max2)
        @test !isempty(fva_min2)
    end
end
