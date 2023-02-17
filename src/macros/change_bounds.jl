"""
    @_change_bounds_fn ModelType IdxType [plural] [inplace] begin ... end

A helper for creating simple bounds-changing function similar to
[`change_bounds`](@ref).
"""
macro _change_bounds_fn(model_type, idx_type, args...)
    body = last(args)
    typeof(body) == Expr || throw(DomainError(body, "missing function body"))
    plural = :plural in args
    plural_s = plural ? "s" : ""
    inplace = :inplace in args
    fname = Symbol(:change_bound, plural_s, inplace ? "!" : "")
    idx_var = Symbol(
        :rxn,
        idx_type == :Int ? "_idx" :
        idx_type == :String ? "_id" :
        throw(DomainError(idx_type, "unsupported index type for change_bound macro")),
        plural_s,
    )
    lower_bound_s = Symbol(:lower_bound, plural_s)
    upper_bound_s = Symbol(:upper_bound, plural_s)

    example_idx =
        plural ? (idx_type == :Int ? [123, 234] : ["ReactionA", "ReactionC"]) :
        (idx_type == :Int ? 123 : "\"ReactionB\"") #= unquoting is hard =#
    example_val = plural ? [4.2, 100.1] : 42.3
    missing_default = plural ? :((nothing for _ in $idx_var)) : nothing

    bound_type = Float64
    if plural
        idx_type = AbstractVector{eval(idx_type)}
        bound_type = AbstractVector{bound_type}
    end

    docstring = """
        $fname(
            model::$model_type,
            $idx_var::$idx_type;
            lower_bound$(plural_s) =$missing_default,
            upper_bound$(plural_s) =$missing_default,
        )

    Change the specified reaction flux bound$(plural_s) in the model
    $(inplace ? "in-place" : "and return the modified model").

    # Example
    ```
    $(inplace ? "new_model = " : "")$fname(model, $example_idx, lower_bound =$(-0.5 .* example_val), upper_bound =$example_val)
    ```
    """

    esc(
        Expr(
            :macrocall,
            Symbol("@doc"),
            __source__,
            docstring,
            :(
                $fname(
                    model::$model_type,
                    $idx_var::$idx_type;
                    $lower_bound_s = $missing_default,
                    $upper_bound_s = $missing_default,
                ) = $body
            ),
        ),
    )
end
