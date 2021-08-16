macro _remove_fn(objname, model_type, idx_type, args...)
    body = last(args)
    typeof(body) == Expr || throw(DomainError(body, "missing function body"))
    plural = :plural in args
    plural_s = plural ? "s" : ""
    inplace = :inplace in args
    fname = Symbol(:remove_, objname, plural_s, inplace ? "!" : "")
    idx_var = Symbol(
        objname,
        idx_type == :Int ? "_idx" :
        idx_type == :String ? "_id" :
        throw(DomainError(idx_type, "unsupported index type for _remove_fn macro")),
        plural_s,
    )

    if plural
        idx_type = AbstractVector{eval(idx_type)}
    end

    docstring = """
        $fname(model::$model_type, $idx_var::$idx_type)

    Remove $objname$plural_s from the model of type `$model_type`
    $(inplace ? "in-place" : "and return the modified model").
    """

    return :(@doc $docstring $fname(model::$model_type, $idx_var::$idx_type) = $body)
end
