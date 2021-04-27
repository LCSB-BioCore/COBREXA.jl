"""
Pretty printing of `Gene`.
"""
function Base.show(io::IO, ::MIME"text/plain", g::Gene)
    for fname in fieldnames(Gene)
        _print_with_colors(io, "Gene.$(string(fname)): ", getfield(g, fname))
    end
end
