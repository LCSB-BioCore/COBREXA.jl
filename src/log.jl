
"""
    module Log

Logging helpers for COBREXA.

# Exports
$(EXPORTS)
"""
module Log
using ..ModuleTools
@dse

"""
    module Internal

Internal helpers for logging.

# Exports
$(EXPORTS)
"""
module Internal
using ..ModuleTools
@dse
"""
$(TYPEDSIGNATURES)

This creates a group of functions that allow masking out topic-related logging
actions. A call that goes as follows:

    @make_logging_tag XYZ

creates the following tools:

- global variable `XYZ_log_enabled` defaulted to false
- function `log_XYZ` that can be called to turn the logging on/off
- a masking macro `@XYZ_log` that can be prepended to commands that should
  only happen if the logging of tag XYZ is enabled.

The masking macro is then used as follows:

    @XYZ_log @info "This is the extra verbose information you wanted!" a b c

The user can direct logging with these:

    log_XYZ()
    log_XYZ(false)

`doc` should be a name of the stuff that is being printed if the corresponding
log_XYZ() is enabled -- it is used to create a friendly documentation for the
logging switch. In this case it could say `"X, Y and Z-related messages"`.
"""
macro make_logging_tag(sym::Symbol, doc::String)
    enable_flag = Symbol(sym, :_log_enabled)
    enable_fun = Symbol(:log_, sym)
    log_macro = Symbol(sym, :_log)
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
                $enable_flag ? x : nothing
            end
        end
    ))
end

@make_logging_tag models "model-related messages"
@make_logging_tag io "messages and warnings from model input/output"
@make_logging_tag perf "performance-related tracing information"

@export_locals
end #Internal

import .Internal: log_models, log_io, log_perf

#TODO can this be exported automatically?
export log_models
export log_io
export log_perf

@export_locals
end
