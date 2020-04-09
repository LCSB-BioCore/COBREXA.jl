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

test_simpleLP2() = LinearModel(zeros(2, 2),
                              [0., 0.],
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

test_coupledLP() = LinearModel(sprand(4000, 3000, 0.5),
                               sprand(4000, 0.5),
                               sprand(2000, 3000, 0.5),
                               sprand(2000, 0.5),
                               sprand(2000, 0.5),
                               sprand(3000, 0.5),
                               sprand(3000, 0.5),
                               sprand(3000, 0.5),
                               ["r$x" for x in 1:3000],
                               ["m$x" for x in 1:4000])