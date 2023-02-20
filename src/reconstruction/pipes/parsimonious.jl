"""
$(TYPEDSIGNATURES)

Pipe-able version of the [`ParsimoniousModel`](@ref) wrapper that minimizes the
solution of variables specified by `semantics`.
"""
with_parsimonious_solution(semantics::Vector{Symbol}) =
    model -> ParsimoniousModel(model, semantics)

"""
$(TYPEDSIGNATURES)

Pipe-able version of the [`ParsimoniousModel`](@ref) wrapper that minimizes the
solution of variables specified by `semantic`.
"""
with_parsimonious_solution(semantic::Symbol) = with_parsimonious_solution([semantic])

"""
$(TYPEDSIGNATURES)

Pipe-able version of the [`ParsimoniousModel`](@ref) wrapper that minimizes the
solution of variables specified by `var_ids`.
"""
with_parsimonious_solution(var_ids::Vector{String}) =
    model -> ParsimoniousModel(model, var_ids)
