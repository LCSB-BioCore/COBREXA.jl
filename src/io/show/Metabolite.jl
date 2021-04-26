
"""
Pretty printing of metabolite::Metabolite.
"""
function Base.show(io::IO, ::MIME"text/plain", m::Metabolite)
    _print_color(io, "Metabolite ID: ", m.id)
    _print_color(io, "Name: ", m.name)
    _print_color(io, "Formula: ", m.formula)
    _print_color(io, "Charge: ", string(m.charge))
    _print_color(io, "Compartment: ", string(m.compartment))
    _print_color(io, "Notes: ", m.notes)
    _print_color(io, "Annotation: ", m.annotation)
    _print_color(io, "Fields: ", join([string(x) for x in fieldnames(Metabolite)], ", "))
end

"""
Pretty printing of metabolites::Vector{Metabolite}.
"""
function Base.show(io::IO, ::MIME"text/plain", ms::Vector{Metabolite})
    _print_color(io, "Metabolite vector of length: : ", string(length(ms)))
    _print_color(
        io,
        "Each metabolite has fields: ",
        join([string(x) for x in fieldnames(Metabolite)], ", "),
    )
end
