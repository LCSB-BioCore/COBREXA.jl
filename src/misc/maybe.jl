
"""
$(TYPEDSIGNATURES)

Helper for getting stuff from dictionaries where keys may be easily missing.
"""
maybeget(_::Nothing, _...) = nothing
maybeget(x, k, ks...) = haskey(x, k) ? maybeget(x[k], ks...) : nothing
maybeget(x) = x

"""
$(TYPEDSIGNATURES)

Helper for applying functions to stuff that might be `nothing`.
"""
maybemap(f, _::Nothing) = nothing
maybemap(f, x) = f(x)
