function Base.show(io::IO, ::MIME"text/plain", g::Gene)
    for fname in fieldnames(Gene)
        _pretty_print_keyvals(io, "Gene.$(string(fname)): ", getfield(g, fname))
    end
end
