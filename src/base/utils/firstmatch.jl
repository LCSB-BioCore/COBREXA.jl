
"""
    _firstmatch(pred, collection, ret=identity)

A helper to linearly find elements in collections by a predicate. Returns
`nothing` if the predicate did not match any element.
"""
function _firstmatch(pred, collection, ret = identity)
    for x in collection
        if pred(x)
            return ret(x)
        end
    end
    return nothing
end
