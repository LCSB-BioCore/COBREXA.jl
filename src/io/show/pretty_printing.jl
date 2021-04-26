"""
    _pretty_print(def::String, payload::String; def_color=:blue, payload_color=:magenta)

Function to faciliate pretty printing of COBREXA structs.
"""
function _print_color(
    io,
    def::String,
    payload::String;
    def_color = :blue,
    payload_color = :magenta,
)
    print(io, Crayon(foreground = def_color), def)
    if isempty(payload)
        println(io, Crayon(foreground = payload_color), "---")
    else
        println(io, Crayon(foreground = payload_color), payload)
    end
end

function _print_color(io, def::String, payload; def_color = :blue, payload_color = :magenta)
    print(io, Crayon(foreground = def_color), def)
    if isempty(payload)
        println(io, Crayon(foreground = payload_color), "---")
    else
        println(io, "")
        for (k, v) in payload
            vv = typeof(v) == String ? v : join(v, ", ") # probably unnecessary once #92 is implemented
            println(io, Crayon(foreground = payload_color), "\t", k, ": ", vv)
        end
    end
end
