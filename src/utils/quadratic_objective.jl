
"""
$(TYPEDSIGNATURES)

Produce a matrix for [`objective`](@ref) so that the model minimizes the
squared distance from the `center` variable assignment. Length of `center`
should be the same as number of variables in the model.
"""
negative_squared_distance_objective(center::Vector{Float64}) =
    [spdiagm(fill(-0.5, length(center))) center]
