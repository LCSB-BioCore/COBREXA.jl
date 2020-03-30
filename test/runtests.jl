using Test
using ***REMOVED***

test_cobraLP() = CobraLP(zeros(4, 3),
                         zeros(4),
                         ones(3),
                         ones(3),
                         ones(3),
                         ["r1"],
                         ["m1"])


@testset "CobraLP type" begin
    cp = test_cobraLP()
    @test cp isa CobraLP
end

@testset "Add reaction" begin
    cp = test_cobraLP()
    cp = addReaction(cp, 2. * ones(3), 1.)
    @test size(cp.S) == (5, 3)
    cp = addReactions(cp, 2. * ones(1, 3), ones(1))
    @test size(cp.S) == (6, 3)
    cp = addReactions(cp, 2. * ones(10, 3), ones(10))
    @test size(cp.S) == (16, 3)
end
