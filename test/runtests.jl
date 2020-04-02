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
                        ["r1"],
                        ["m1"])

test_simpleLP() = LinearModel([ 1.0 1.0
                               -1.0 1.0],
                              [3., 1.],
                              [-0.25, 1.],
                              -ones(2),
                              2.0 * ones(2),
                              ["r1"],
                              ["m1"])

test_sparseLP() = LinearModel(sprand(4000, 3000, 0.5),
                              sprand(4000, 0.5),
                              sprand(3000, 0.5),
                              sprand(3000, 0.5),
                              sprand(3000, 0.5),
                              ["r1"],
                              ["m1"])

@testset "LinearModel type" begin
    cp = test_LP()
    @test cp isa LinearModel
end

@testset "Add reactions" begin
    cp = test_LP()
    @test size(cp.S) == (4, 3)
    cp = addReactions(cp, 2. * ones(3), 1.)
    @test size(cp.S) == (5, 3)
    cp = addReactions(cp, 2. * ones(1, 3), ones(1))
    @test size(cp.S) == (6, 3)
    cp = addReactions(cp, 2. * ones(10, 3), ones(10))
    @test size(cp.S) == (16, 3)

    cp = test_sparseLP()
    @test size(cp.S) == (4000, 3000)
    cp = addReactions(cp, 2. * sprand(3000, 0.5), 1.)
    @test size(cp.S) == (4001, 3000)
    cp = addReactions(cp, 2. * sprand(1, 3000, 0.5), sprand(1, 0.5))
    @test size(cp.S) == (4002, 3000)
    cp = addReactions(cp, 2. * sprand(1000, 3000, 0.5), sprand(1000, 0.5))
    @test size(cp.S) == (5002, 3000)

    cp = test_sparseLP()
    @test size(cp.S) == (4000, 3000)
    cp = addReactions(cp, 2. * ones(3000), 1.)
    @test size(cp.S) == (4001, 3000)
    cp = addReactions(cp, 2. * ones(1, 3000), ones(1))
    @test size(cp.S) == (4002, 3000)
    cp = addReactions(cp, 2. * ones(1000, 3000), ones(1000))
    @test size(cp.S) == (5002, 3000)
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
    @test_throws ErrorException loadModel("agora-model.mat", "badmodel")
end
