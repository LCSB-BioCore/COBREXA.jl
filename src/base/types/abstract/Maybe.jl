
"""
    Maybe{T} = Union{Nothing, T}

A nice name for "nullable" type.
"""
const Maybe{T} = Union{Nothing,T}

"""
    default(d::T, x::Maybe{T})::T where {T}

Fold the `Maybe{T}` down to `T` by defaulting.
"""
function default(d::T, x::Maybe{T})::T where {T}
    isnothing(x) ? d : x
end

function maybemap(f, x::Maybe)::Maybe
    isnothing(x) ? nothing : f(x)
end
