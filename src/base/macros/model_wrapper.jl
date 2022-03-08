"""
    @_inherit_model_methods

Generates trivial accessor functions listed in `fns` for a model that is
wrapped in type `mtype` as field `member`.
"""
macro _inherit_model_methods(mtype::Symbol, arglist, member::Symbol, fwdlist, fns...)
    Expr(
        :block,
        (
            begin
                header = Expr(:call, fn, :(model::$mtype), arglist.args...)
                call = Expr(:call, fn, :(model.$member), fwdlist.args...)
                esc(
                    Expr(
                        :macrocall,
                        Symbol("@doc"),
                        __source__,
                        """
                            $header

                        Evaluates [`$fn`](@ref) on the model contained in $mtype.
                        """,
                        Expr(:(=), header, Expr(:block, __source__, call)),
                    ),
                )
            end for fn in fns
        )...,
    )
end
