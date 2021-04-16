test_LP() = CoreModel(
    zeros(4, 3),
    zeros(4),
    ones(3),
    ones(3),
    ones(3),
    ["r$x" for x = 1:3],
    ["m$x" for x = 1:4],
)

test_simpleLP() = CoreModel(
    [
        1.0 1.0
        -1.0 1.0
    ],
    [3.0, 1.0],
    [-0.25, 1.0],
    -ones(2),
    2.0 * ones(2),
    ["r$x" for x = 1:2],
    ["m$x" for x = 1:2],
)

test_simpleLP2() = CoreModel(
    zeros(2, 2),
    [0.0, 0.0],
    [-0.25, 1.0],
    -ones(2),
    2.0 * ones(2),
    ["r$x" for x = 1:2],
    ["m$x" for x = 1:2],
)

test_sparseLP() = CoreModel(
    sprand(4000, 3000, 0.5),
    sprand(4000, 0.5),
    sprand(3000, 0.5),
    sprand(3000, 0.5),
    sprand(3000, 0.5),
    ["r$x" for x = 1:3000],
    ["m$x" for x = 1:4000],
)

test_coupledLP() = CoreModelCoupled(
    CoreModel(
        sprand(4000, 3000, 0.5),
        sprand(4000, 0.5),
        sprand(3000, 0.5),
        sprand(3000, 0.5),
        sprand(3000, 0.5),
        ["r$x" for x = 1:3000],
        ["m$x" for x = 1:4000],
    ),
    sprand(2000, 3000, 0.5),
    sprand(2000, 0.5),
    sprand(2000, 0.5),
)
