"""
    _serialized_change_unwrap(fn::Symbol)

Creates a simple wrapper structure that calls a function transparently on the
internal precached model. Internal type is returned (because this would break
the consistency of serialization).
"""
macro _serialized_change_unwrap(fn::Symbol)
    docstring = """
        $fn(model::Serialized, ...)

    Calls [$fn](@ref) of the internal serialized model type.
    Returns the modified un-serialized model.
    """
    :(@doc $docstring $fn(model::Serialized, args...; kwargs...) =
        $fn(unwrap_serialized(model), args...; kwargs...))
end
