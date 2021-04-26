"""
Pretty printing of `Gene`.
"""
function Base.show(io::IO, ::MIME"text/plain", g::Gene)
    for fname in fieldnames(Gene)
        _print_color(io, "Gene.$(string(fname)): ", getfield(g, fname))
    end
end

"""
Pretty printing of gene reaction rules in type `Vector{Vector{Gene}}`.
"""
function Base.show(io::IO, ::MIME"text/plain", grr::Vector{Vector{Gene}})
    grr_strings = String[]
    for gr in grr
        push!(grr_strings, "(" * join([g.id for g in gr], " and ") * ")")
    end
    _print_color(io, "Gene reaction rule: ", join(grr_strings, " or "))
end
