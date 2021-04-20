
"""
Pretty printing of metabolite::Metabolite.
"""
function Base.show(io::IO, ::MIME"text/plain", m::Metabolite)
    _pretty_print(io, "Metabolite ID: ", m.id)
    _pretty_print(io, "Name: ", m.name)
    _pretty_print(io, "Formula: ", m.formula)
    _pretty_print(io, "Charge: ", string(m.charge))
    _pretty_print(io, "Compartment: ", string(m.compartment))
    _pretty_print(io, "Notes: ", m.notes)
    _pretty_print(io, "Annotation: ", m.annotation)
    _pretty_print(io, "Fields: ", join([string(x) for x in fieldnames(Metabolite)], ", "))
end

"""
Pretty printing of metabolites::Vector{Metabolite}.
"""
function Base.show(io::IO, ::MIME"text/plain", ms::Vector{Metabolite})
    _pretty_print(io, "Metabolite vector of length: : ", string(length(ms)))
    _pretty_print(
        io,
        "Each metabolite has fields: ",
        join([string(x) for x in fieldnames(Metabolite)], ", "),
    )
end
