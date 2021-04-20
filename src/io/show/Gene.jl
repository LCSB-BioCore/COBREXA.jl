"""
Pretty printing of `Gene`.
"""
function Base.show(io::IO, ::MIME"text/plain", g::Gene)
    _pretty_print(io, "Gene ID: ", g.id)
    _pretty_print(io, "Name: ", g.name)
    _pretty_print(io, "Notes: ", g.notes)
    _pretty_print(io, "Annotation: ", g.annotation)
    _pretty_print(io, "Fields: ", join([string(x) for x in fieldnames(Gene)],", "))
end

"""
Pretty printing of `Vector{Gene}`.
"""
function Base.show(io::IO, ::MIME"text/plain", gs::Vector{Gene})
    _pretty_print(io, "Gene vector of length: ", string(length(gs)))
    _pretty_print(io, "Each gene has fields: ", join([string(x) for x in fieldnames(Gene)], ", "))
end

"""
Pretty printing of gene reaction rules in type `Vector{Vector{Gene}}`.
"""
function Base.show(io::IO, ::MIME"text/plain", grr::Vector{Vector{Gene}})
    grr_strings = String[]
    for gr in grr
        push!(grr_strings, "(" * join([g.id for g in gr], " and ") * ")")
    end
    _pretty_print(io, "Gene reaction rule: ", join(grr_strings, " or "))
end
