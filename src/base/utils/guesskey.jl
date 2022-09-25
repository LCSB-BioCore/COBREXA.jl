"""
$(TYPEDSIGNATURES)

Unfortunately, many model types that contain dictionares do not have
standardized field names, so we need to try a few possibilities and guess the
best one. The keys used to look for valid field names should be ideally
specified as constants in `src/base/constants.jl`.
"""
function _guesskey(avail, possibilities)
    x = intersect(possibilities, avail)

    if isempty(x)
        @debug "could not find any of keys: $possibilities"
        return nothing
    end

    if length(x) > 1
        @debug "Possible ambiguity between keys: $x"
    end
    return x[1]
end

"""
$(TYPEDSIGNATURES)

Return `fail` if key in `keys` is not in `collection`, otherwise
return `collection[key]`. Useful if may different keys need to be
tried due to non-standardized model formats.
"""
function _gets(collection, fail, keys)
    for key in keys
        haskey(collection, key) && return collection[key]
    end
    return fail
end