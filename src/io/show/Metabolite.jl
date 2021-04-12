
"""
Pretty printing of metabolite::Metabolite.
"""
function Base.show(io::IO, ::MIME"text/plain", m::Metabolite)
    println(
        io,
        "Metabolite ID: ",
        m.id,
        "\n",
        "Metabolite name: ",
        m.name,
        "\n",
        "Formula: ",
        m.formula,
        "\n",
        "Charge: ",
        m.charge,
    )
end

"""
Pretty printing of metabolites::Vector{Metabolite}.
"""
function Base.show(io::IO, ::MIME"text/plain", ms::Vector{Metabolite})
    println(io, "Metabolite set of length: ", length(ms))
end
