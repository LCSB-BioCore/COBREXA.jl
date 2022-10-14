
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

Apply a function to `x` only if it is not `nothing`.
"""
function _maybemap(f, x::Maybe)::Maybe
    isnothing(x) ? nothing : f(x)
end
