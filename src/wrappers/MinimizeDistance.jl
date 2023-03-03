
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
    model -> MinimizeSolutionDistance(center, model)

"""
$(TYPEDEF)

A wrapper that sets an objective that minimizes Euclidean distance from a given
point in a semantic.
"""
struct MinimizeSemanticDistance <: AbstractModelWrapper
    semantic::Symbol
    center::Vector{Float64}
    inner::AbstractMetabolicModel
end

Accessors.unwrap_model(m::MinimizeSemanticDistance) = m.inner

Accessors.objective(m::MinimizeSemanticDistance) =
    let
        Sem = variable_semantic_mtx(m.semantic, m.inner)
        Sem' *
        [spdiagm(fill(-0.5, size(Sem, 2))) m.center] *
        [Sem zeros(size(Sem, 1)); zeros(size(Sem, 2))' 1.0]
    end

"""
$(TYPEDSIGNATURES)

Set a quadratic objective that minimizes the solution distance from a selected
point in a space defined by a given semantic. Use
[`minimize_projected_distance`](@ref) for more fine-grained and weighted
variant, and [`minimize_solution_distance`](@ref) for working directly upon
variables. Internally powered by [`MinimizeSemanticDistance`](@ref).
"""
minimize_semantic_distance(semantic::Symbol, center::Vector{Float64}) =
    model -> MinimizeSolutionDistance(semantic, center, model)

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
    model -> MinimizeProjectedDistance(proj, center, model)
