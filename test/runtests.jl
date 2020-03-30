using Test
using ***REMOVED***

test_cobraLP() = CobraLP([1.0 0.],
                         [0.],
                         [0.],
                         [0.],
                         [0.],
                         ["r1"],
                         ["m1"])


@testset "CobraLP type" begin
    cp = test_cobraLP()
    @test cp isa CobraLP
end
