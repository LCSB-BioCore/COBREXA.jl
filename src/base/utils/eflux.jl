
"""
$(TYPEDSIGNATURES)

Create E-Flux-like bounds for a reaction that requires gene products given by
`grr`, with expression "choke" data given by relative probabilities (between 0
and 1) in `relative_expression`.
"""
function expression_probabilistic_bounds(
    relative_expression::Dict{String,Float64},
    grr::Maybe{Vector{Vector{String}}},
)::Float64
    isnothing(grr) && return 1.0
    lup(g) = get(relative_expression, g, 1.0)
    1 - prod(1 - prod(lup.(gr)) for gr in grr)
end
