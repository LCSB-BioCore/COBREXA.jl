
"""
    _firstmatch(pred, collection, ret=identity)

A helper to linearly find elements in collections by a predicate. Returns
`nothing` if the predicate did not match any element.
"""
function _firstmatch(pred, collection, ret = identity)
    for i in collection
        if pred(i)
            return ret(i)
        end
    end
    return nothing
end
