"""
    _pretty_print(io, def::String, payload::String; def_color=:cyan, empty_color=:dark_gray, payload_color=:default)

Function to faciliate pretty printing of COBREXA structs.
"""
function _print_color(
    io,
    def::String,
    payload::String;
    def_color = _constants.color.key,
    empty_color = _constants.color.empty,
    payload_color = _constants.color.payload,
)
    print(io, Crayon(foreground = def_color), def)
    if isempty(payload)
        println(io, Crayon(foreground = empty_color), "---")
    else
        println(io, Crayon(foreground = payload_color), payload)
    end
end

function _print_color(
    io,
    def::String,
    payload::Dict;
    def_color = _constants.color.key,
    empty_color = _constants.color.empty,
    payload_color = _constants.color.payload,
)

    print(io, Crayon(foreground = def_color), def)
    if isempty(payload)
        println(io, Crayon(foreground = empty_color), "---")
    else
        println(io, "")
        for (k, v) in payload
            vv = typeof(v) == String ? v : join(v, ", ") # probably unnecessary once #92 is implemented
            println(io, Crayon(foreground = payload_color), "\t", k, ": ", vv)
        end
    end
end
