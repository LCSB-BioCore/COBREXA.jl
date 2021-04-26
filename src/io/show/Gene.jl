"""
Pretty printing of `Gene`.
"""
function Base.show(io::IO, ::MIME"text/plain", g::Gene)
    _print_color(io, "Gene ID: ", g.id)
    _print_color(io, "Name: ", g.name)
    _print_color(io, "Notes: ", g.notes)
    _print_color(io, "Annotation: ", g.annotation)
    _print_color(io, "Fields: ", join([string(x) for x in fieldnames(Gene)], ", "))
end

"""
Pretty printing of `Vector{Gene}`.
"""
function Base.show(io::IO, ::MIME"text/plain", gs::Vector{Gene})
    _print_color(io, "Gene vector of length: ", string(length(gs)))
    _print_color(
        io,
        "Each gene has fields: ",
        join([string(x) for x in fieldnames(Gene)], ", "),
    )
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
