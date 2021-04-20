"""
    info_msg(str, atlevel=10)

Only print this message `str` if `env[\"log-level\"]` ≥`atlevel`.
"""
function info_msg(str, atlevel = 10)
    if COBREXA.env["log-level"] >= atlevel
        println(Crayon(foreground = :cyan), str)
    end
end

"""
    warn_msg(str, atlevel=10)

Only print this message `str` if `env[\"log-level\"]` ≥`atlevel`.
"""
function warn_msg(str, atlevel = 10)
    if COBREXA.env["log-level"] >= atlevel
        println(Crayon(foreground = :red), str)
    end
end
