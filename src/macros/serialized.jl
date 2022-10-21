"""
    @_serialized_change_unwrap function

Creates a simple wrapper structure that calls the `function` transparently on
the internal precached model. The internal type is returned (otherwise this
would break the consistency of serialization).
"""
macro _serialized_change_unwrap(fn::Symbol)
    docstring = """
        $fn(model::Serialized, ...)

    Calls [`$fn`](@ref) of the internal serialized model type.
    Returns the modified un-serialized model.
    """
    esc(
        Expr(
            :macrocall,
            Symbol("@doc"),
            __source__,
            docstring,
            :(
                $fn(model::Serialized, args...; kwargs...) =
                    $fn(unwrap_serialized(model), args...; kwargs...)
            ),
        ),
    )
end
