"""
Pretty printing of `Gene`.
"""
function Base.show(io::IO, ::MIME"text/plain", g::Gene)
    for fname in fieldnames(Gene)
        _print_with_colors(io, "Gene.$(string(fname)): ", getfield(g, fname))
    end
end

"""
Pretty printing of gene reaction rules in type `Vector{Vector{Gene}}`.
"""
function Base.show(io::IO, ::MIME"text/plain", grr::Vector{Vector{Gene}})
    _print_with_colors(
        io,
        "Gene reaction rule: ",
        _unparse_grr(grr),
    )
end
