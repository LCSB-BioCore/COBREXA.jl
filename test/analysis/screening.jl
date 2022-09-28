
@testset "Screening functions" begin
    m = test_toyModel()

    # nothing to analyze
    @test_throws DomainError screen(m; analysis = identity)

    # array dimensionalities must match
    @test_throws DomainError screen(
        m;
        analysis = identity,
        variants = [[], []],
        args = [() ()],
    )

    # sizes must match
    @test_throws DomainError screen(
        m;
        analysis = identity,
        variants = [[], []],
        args = [()],
    )

    # argument handling
    @test screen(m, analysis = identity, variants = [[]]) == [m]
    @test screen(m, analysis = identity, args = [()]) == [m]
    @test screen(m, analysis = (a, b) -> b, args = [(1,), (2,)]) == [1, 2]

    # test modifying some reactions
    quad_rxn(i) = (m::CoreModel) -> begin
        mm = copy(m)
        mm.S = copy(m.S)
        mm.S[:, i] .^= 2
        return mm
    end

    @test screen_variants(
        m,
        [[quad_rxn(i)] for i = 1:3],
        m -> Analysis.flux_balance_analysis_vec(m, Tulip.Optimizer);
        workers = W,
    ) == [
        [250.0, -250.0, -1000.0, 250.0, 1000.0, 250.0, 250.0],
        [500.0, 500.0, 1000.0, 500.0, -1000.0, 500.0, 500.0],
        [500.0, 500.0, 1000.0, -500.0, 1000.0, 500.0, 500.0],
    ]

    # test solver modifications
    @test screen(
        m;
        analysis = (m, sense) -> Analysis.flux_balance_analysis_vec(
            m,
            Tulip.Optimizer;
            modifications = [change_sense(sense)],
        ),
        args = [(MIN_SENSE,), (MAX_SENSE,)],
    ) == [
        [-500.0, -500.0, -1000.0, 500.0, 1000.0, -500.0, -500.0],
        [500.0, 500.0, 1000.0, -500.0, -1000.0, 500.0, 500.0],
    ]
end
