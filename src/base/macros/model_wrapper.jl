
"""
$(TYPEDSIGNATURES)

A helper backend for [`@_inherit_model_methods`](@ref) and
[`@_inherit_model_methods_fn`](@ref).
"""
function _inherit_model_methods_impl(
    source,
    mtype::Symbol,
    arglist,
    access,
    fwdlist,
    fns...,
)
    Expr(
        :block,
        (
            begin
                header = Expr(:call, fn, :(model::$mtype), arglist.args...)
                call = Expr(:call, fn, access(:model), fwdlist.args...)
                esc(
                    Expr(
                        :macrocall,
                        Symbol("@doc"),
                        source,
                        """
                            $header

                        Evaluates [`$fn`](@ref) on the model contained in $mtype.
                        """,
                        Expr(:(=), header, Expr(:block, source, call)),
                    ),
                )
            end for fn in fns
        )...,
    )
end

"""
    @_inherit_model_methods

Generates trivial accessor functions listed in `fns` for a model that is
wrapped in type `mtype` as field `member`.
"""
macro _inherit_model_methods(mtype::Symbol, arglist, member::Symbol, fwdlist, fns...)
    _inherit_model_methods_impl(
        __source__,
        mtype,
        arglist,
        sym -> :($sym.$member),
        fwdlist,
        fns...,
    )
end

"""
    @_inherit_model_methods_fn

A more generic version of [`@_inherit_model_methods`](@ref) that accesses the
"inner" model using an accessor function name.
"""
macro _inherit_model_methods_fn(mtype::Symbol, arglist, accessor, fwdlist, fns...)
    _inherit_model_methods_impl(
        __source__,
        mtype,
        arglist,
        sym -> :($accessor($sym)),
        fwdlist,
        fns...,
    )
end
