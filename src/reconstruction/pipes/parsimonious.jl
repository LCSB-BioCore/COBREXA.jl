"""
$(TYPEDSIGNATURES)

Pipe-able version of the [`ParsimoniousModel`](@ref) wrapper.
"""
with_parsimonious_solution(semantics::Vector{Symbol}) =
    model -> ParsimoniousModel(model, semantics)

with_parsimonious_solution(semantic::Symbol) = with_parsimonious_solution([semantic])
