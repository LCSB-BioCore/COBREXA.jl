
"""
$(TYPEDSIGNATURES)

Fold the `Maybe{T}` down to `T` by defaulting.
"""
function default(d::T, x::Maybe{T})::T where {T}
    isnothing(x) ? d : x
end

"""
$(TYPEDSIGNATURES)

Apply a function to `x` only if it is not `nothing`.
"""
function maybemap(f, x::Maybe)::Maybe
    isnothing(x) ? nothing : f(x)
end
