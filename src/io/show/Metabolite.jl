function Base.show(io::Base.IO, ::MIME"text/plain", m::Metabolite)
    for fname in fieldnames(Metabolite)
        if fname == :charge
            c = isnothing(getfield(m, fname)) ? nothing : string(getfield(m, fname))
            _pretty_print_keyvals(io, "Metabolite.$(string(fname)): ", c)
        else
            _pretty_print_keyvals(io, "Metabolite.$(string(fname)): ", getfield(m, fname))
        end
    end
end
