using Test
using SparseArrays
using ***REMOVED***

test_cobraLP() = CobraLP(zeros(4, 3),
                         zeros(4),
                         ones(3),
                         ones(3),
                         ones(3),
                         ["r1"],
                         ["m1"])

test_sparseCobraLP() = CobraLP(sprand(4000, 3000, 0.5),
                               sprand(4000, 0.5),
                               sprand(3000, 0.5),
                               sprand(3000, 0.5),
                               sprand(3000, 0.5),
                               ["r1"],
                               ["m1"])

@testset "CobraLP type" begin
    cp = test_cobraLP()
    @test cp isa CobraLP
end

@testset "Add reactions" begin
    cp = test_cobraLP()
    @test size(cp.S) == (4, 3)
    cp = addReaction(cp, 2. * ones(3), 1.)
    @test size(cp.S) == (5, 3)
    cp = addReactions(cp, 2. * ones(1, 3), ones(1))
    @test size(cp.S) == (6, 3)
    cp = addReactions(cp, 2. * ones(10, 3), ones(10))
    @test size(cp.S) == (16, 3)

    cp = test_sparseCobraLP()
    @test size(cp.S) == (4000, 3000)
    cp = addReaction(cp, 2. * sprand(3000, 0.5), 1.)
    @test size(cp.S) == (4001, 3000)
    cp = addReactions(cp, 2. * sprand(1, 3000, 0.5), sprand(1, 0.5))
    @test size(cp.S) == (4002, 3000)
    cp = addReactions(cp, 2. * sprand(1000, 3000, 0.5), sprand(1000, 0.5))
    @test size(cp.S) == (5002, 3000)

    cp = test_sparseCobraLP()
    @test size(cp.S) == (4000, 3000)
    cp = addReaction(cp, 2. * ones(3000), 1.)
    @test size(cp.S) == (4001, 3000)
    cp = addReactions(cp, 2. * ones(1, 3000), ones(1))
    @test size(cp.S) == (4002, 3000)
    cp = addReactions(cp, 2. * ones(1000, 3000), ones(1000))
    @test size(cp.S) == (5002, 3000)
end
