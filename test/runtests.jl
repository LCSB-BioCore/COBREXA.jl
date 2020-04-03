using Test
using SparseArrays
using ***REMOVED***
using JuMP
using MAT

test_LP() = LinearModel(zeros(4, 3),
                        zeros(4),
                        ones(3),
                        ones(3),
                        ones(3),
                        ["r$x" for x in 1:3],
                        ["m$x" for x in 1:4])

test_simpleLP() = LinearModel([ 1.0 1.0
                               -1.0 1.0],
                              [3., 1.],
                              [-0.25, 1.],
                              -ones(2),
                              2.0 * ones(2),
                              ["r$x" for x in 1:2],
                              ["m$x" for x in 1:2])

test_sparseLP() = LinearModel(sprand(4000, 3000, 0.5),
                              sprand(4000, 0.5),
                              sprand(3000, 0.5),
                              sprand(3000, 0.5),
                              sprand(3000, 0.5),
                              ["r$x" for x in 1:3000],
                              ["m$x" for x in 1:4000])

@testset "LinearModel type" begin
    cp = test_LP()
    @test cp isa LinearModel
end

@testset "Add reactions" begin
    cp = test_LP()
    @test size(cp.S) == (4, 3)
    cp = addReactions(cp, 2. * ones(4), 2., -1., 1.)
    @test size(cp.S) == (4, 4)
    cp = addReactions(cp, 2. * ones(4, 1), 2 .* ones(1), -ones(1), ones(1))
    @test size(cp.S) == (4, 5)
    cp = addReactions(cp, 2. * ones(4, 10), 2 .* ones(10), -ones(10), ones(10))
    @test size(cp.S) == (4, 15)

    cp = test_sparseLP()
    @test size(cp.S) == (4000, 3000)
    cp = addReactions(cp, 2. * sprand(4000, 0.5), 2., -1., 1.)
    @test size(cp.S) == (4000, 3001)
    cp = addReactions(cp, 2. * sprand(4000, 1, 0.5), 2 .* sprand(1, 0.5), -sprand(1, 0.5), sprand(1, 0.5))
    @test size(cp.S) == (4000, 3002)
    cp = addReactions(cp, 2. * sprand(4000, 1000, 0.5), 2 .* sprand(1000, 0.5), -sprand(1000, 0.5), sprand(1000, 0.5))
    @test size(cp.S) == (4000, 4002)

    cp = test_sparseLP()
    @test size(cp.S) == (4000, 3000)
    cp = addReactions(cp, 2. * ones(4000), 2., -1., 1.)
    @test size(cp.S) == (4000, 3001)
    cp = addReactions(cp, 2. * ones(4000, 1), 2 .* ones(1), -ones(1), ones(1))
    @test size(cp.S) == (4000, 3002)
    cp = addReactions(cp, 2. * ones(4000, 1000), 2 .* ones(1000), -ones(1000), ones(1000))
    @test size(cp.S) == (4000, 4002)
end


@testset "Solve LP" begin
    cp = test_simpleLP()
    (lp, x) = solveLP(cp)
    @test termination_status(lp) === MOI.OPTIMAL
    sol = JuMP.value.(x)
    @test sol ≈ [1., 2.]
end

@testset "Import MAT" begin
    cp = loadModel("agora-model.mat", "model")
    @test cp isa LinearModel
    @test size(cp.S) == (475, 496)
    @test_throws ErrorException loadModel("agora-model.mat", "badmodel")
end
