
"""
    Maybe{T} = Union{Nothing, T}

A nice name for "nullable" type.
"""
const Maybe{T} = Union{Nothing,T}

"""
$(TYPEDSIGNATURES)

Fold the `Maybe{T}` down to `T` by defaulting.
"""
function _default(d::T, x::Maybe{T})::T where {T}
    isnothing(x) ? d : x
end

"""
$(TYPEDSIGNATURES)

Apply a function to `x` only if it is not `nothing`. Returns `f(x)` when
applied, otherwise returns `default`.
"""
function _maybemap(f, x::Maybe, default = nothing)::Maybe
    isnothing(x) ? default : f(x)
end
