
"""
    macro _make_logging_group(sym::Symbol, doc::String)

This creates a group of functions that allow masking out topic-related logging
actions. A call that goes as follows:

    @_make_logging_tag XYZ

creates the following tools:

- global variable `_XYZ_log_enabled` defaulted to false
- function `log_XYZ` that can be called to turn the logging on/off
- a masking macro `@_XYZ_log` that can be prepended to commands that should
  only happen if the logging of tag XYZ is enabled.

The masking macro is then used as follows:

    @_XYZ_log @info "This is the extra verbose information you wanted!" a b c

The user can direct logging with these:

    log_XYZ()
    log_XYZ(false)

`doc` should be a name of the stuff that is being printed if the corresponding
log_XYZ() is enabled -- it is used to create a friendly documentation for the
logging switch. In this case it could say `"X, Y and Z-related messages"`.
"""
macro _make_logging_tag(sym::Symbol, doc::String)
    enable_flag = Symbol(:_, sym, :_log_enabled)
    enable_fun = Symbol(:log_, sym)
    log_macro = Symbol(:_, sym, :_log)
    # esc() is necessary here because the internal macro processing would
    # otherwise bind the variables incorrectly.
    esc(:(
        begin
            $enable_flag = false

            """
                $(string($enable_fun))(enable::Bool=true)

            Enable (default) or disable (by passing `false`) output of $($doc).
            """
            $enable_fun(x::Bool = true) = begin
                global $enable_flag = x
            end

            macro $log_macro(x)
                $enable_flag ? x : :nothing
            end
        end
    ))
end

@_make_logging_tag models "model-related messages"
@_make_logging_tag io "messages and warnings from model input/output"
@_make_logging_tag perf "performance-related tracing information"
@_make_logging_tag utilities "Utility function related logging information"
log_utilities(true) # switch log on
