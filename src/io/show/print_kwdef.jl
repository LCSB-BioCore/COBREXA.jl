
"""
$(TYPEDSIGNATURES)

Prints out a display-style (multi-line etc.) version of the structure `x`
defined using `Base.@kwdef` to the `io` stream. Falls back to
[`print_kwdef`](@ref) if the IO context has the `:compact` flag set.
"""
function display_kwdef(io::Base.IO, x::T) where {T}
    if get(io, :compact, false)
        print_kwdef(io, x)
    else
        print(io, T)
        println(io, "(")
        for field in fieldnames(T)
            print(io, " ")
            print(io, field)
            print(io, " = ")
            show(IOContext(io, :compact => true), getfield(x, field))
            println(io, ",")
        end
        print(io, ")")
    end
    nothing
end

"""
$(TYPEDSIGNATURES)

Prints out an inline-style (single-line) version of the structure `x` defined
using `Base.@kwdef` to the `io` stream.
"""
function print_kwdef(io::Base.IO, x::T) where {T}
    print(io, T)
    print(io, "(")
    first = true
    for field in fieldnames(T)
        if first
            first = false
        else
            print(io, ", ")
        end
        print(io, field)
        print(io, "=")
        show(IOContext(io, :compact => true), getfield(x, field))
    end
    print(io, ")")
    nothing
end

"""
$(TYPEDSIGNATURES)

Creates overloads of `Base.show` and `Base.print` for a given type.
"""
macro kwdef_printing(t)
    :(
        begin
            Base.show(io::Base.IO, ::MIME"text/plain", x::$t) = display_kwdef(io, x)
            Base.show(io::Base.IO, x::$t) = print_kwdef(io, x)
            Base.print(io::Base.IO, x::$t) = print_kwdef(io, x)
        end
    )
end

@kwdef_printing Reaction
@kwdef_printing Metabolite
@kwdef_printing Gene
@kwdef_printing Isozyme
