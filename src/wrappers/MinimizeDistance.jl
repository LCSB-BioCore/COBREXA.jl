
"""
$(TYPEDEF)

A wrapper for setting an Euclidean-distance-from-a-reference-point-minimizing
objective on solution variables.
"""
struct MinimizeSolutionDistance <: AbstractModelWrapper
    center::Vector{Float64}
    inner::AbstractMetabolicModel
end

Accessors.unwrap_model(m::MinimizeSolutionDistance) = m.inner

Accessors.objective(m::MinimizeSolutionDistance) =
    [spdiagm(fill(-0.5, length(m.center))) m.center]

"""
$(TYPEDSIGNATURES)

Set a quadratic objective that minimizes the solution distance from a selected
point.  Use [`minimize_semantic_distance`](@ref) or
[`minimize_projected_distance`](@ref) for more fine-grained and weighted
variants of the same concept. Internally powered by
[`MinimizeSolutionDistance`](@ref).
"""
minimize_solution_distance(center::Vector{Float64}) =
    model::AbstractMetabolicModel -> MinimizeSolutionDistance(center, model)

"""
$(TYPEDSIGNATURES)

Set an objective that finds a solution of minimal norm. Typically, it is
necessary to also add more constraints to the objective that prevent optimality
of the trivial zero solution.

This can be used to implement [`parsimonious_flux_balance_analysis`](@ref) in a
flexible way that fits into larger model systems.
"""
with_parsimonious_objective() =
    model::AbstractMetabolicModel ->
        MinimizeSolutionDistance(zeros(variable_count(model)), model)

"""
$(TYPEDEF)

A wrapper that sets an objective that minimizes Euclidean distance from a given
point in a semantics.
"""
struct MinimizeSemanticDistance <: AbstractModelWrapper
    semantics::Symbol
    center::Vector{Float64}
    inner::AbstractMetabolicModel
end

Accessors.unwrap_model(m::MinimizeSemanticDistance) = m.inner

function Accessors.objective(m::MinimizeSemanticDistance)
    s = Accessors.Internal.semantics(m.semantics)
    M = s.mapping_matrix(m.inner)

    return M *
           [spdiagm(fill(-0.5, size(M, 2))) m.center] *
           [M' zeros(size(M, 2)); zeros(size(M, 1))' 1.0]
end

"""
$(TYPEDSIGNATURES)

Set a quadratic objective that minimizes the solution distance from a selected
point in a space defined by a given semantics. Use
[`minimize_projected_distance`](@ref) for more fine-grained and weighted
variant, and [`minimize_solution_distance`](@ref) for working directly upon
variables. Internally powered by [`MinimizeSemanticDistance`](@ref).
"""
minimize_semantic_distance(semantics::Symbol, center::Vector{Float64}) =
    model::AbstractMetabolicModel -> MinimizeSolutionDistance(semantics, center, model)

"""
$(TYPEDSIGNATURES)

Set an objective that finds a solution of minimal norm in a given semantics.
This can be used to implement various realistic variants of
[`parsimonious_flux_balance_analysis`](@ref).
"""
with_parsimonious_objective(semantics::Symbol) =
    model::AbstractMetabolicModel -> let
        s = Accessors.Internal.semantics(semantics)
        MinimizeSemanticDistance(semantics, zeros(s.count(model)), model)
    end

"""
$(TYPEDEF)

A wrapper that sets an objective that minimizes Euclidean distance from a given
point in a space defined by a projection matrix.
"""
struct MinimizeProjectedDistance <: AbstractModelWrapper
    proj::SparseMat
    center::Vector{Float64}
    inner::AbstractMetabolicModel
end

Accessors.unwrap_model(m::MinimizeProjectedDistance) = m.inner

Accessors.objective(m::MinimizeProjectedDistance) =
    let
        m.proj' *
        [spdiagm(fill(-0.5, length(m.center))) m.center] *
        [m.proj zeros(size(R, 1)); zeros(size(R, 2))' 1.0]
    end

"""
$(TYPEDSIGNATURES)

Set a quadratic objective that minimizes the solution distance from a selected
point in a space defined by a custom projection matrix. See
[`minimize_solution_distance`](@ref) and [`minimize_semantic_distance`](@ref)
for simpler variants. Internally powered by
[`MinimizeProjectedDistance`](@ref).
"""
minimize_projected_distance(proj::SparseMat, center::Vector{Float64}) =
    model::AbstractMetabolicModel -> MinimizeProjectedDistance(proj, center, model)

"""
$(TYPEDSIGNATURES)

Set a quadratic objective that minimizes the norm in a given projection of
model variables.
"""
with_parsimonious_objective(proj::SparseMat) =
    model::AbstractMetabolicModel ->
        MinimizeProjectedDistance(proj, zeros(size(proj, 1)), model)
