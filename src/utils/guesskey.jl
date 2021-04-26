
"""
    _guesskey(ks, possibilities)

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
