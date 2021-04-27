
"""
Pretty printing of metabolite::Metabolite.
"""
function Base.show(io::IO, ::MIME"text/plain", m::Metabolite)
    for fname in fieldnames(Metabolite)
        if fname == :charge
            _print_with_colors(io, "Metabolite.$(string(fname)): ", string(getfield(m, fname)))
        else
            _print_with_colors(io, "Metabolite.$(string(fname)): ", getfield(m, fname))
        end
    end
end
