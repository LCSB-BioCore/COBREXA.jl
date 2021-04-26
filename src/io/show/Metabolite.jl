
"""
Pretty printing of metabolite::Metabolite.
"""
function Base.show(io::IO, ::MIME"text/plain", m::Metabolite)
    for fname in fieldnames(Metabolite)
        if fname == :charge
            _print_color(io, "Metabolite.$(string(fname)): ", string(getfield(m, fname)))
        else
            _print_color(io, "Metabolite.$(string(fname)): ", getfield(m, fname))
        end
    end
end
