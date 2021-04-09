
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
Pretty printing of metabolites::Array{Metabolite, 1}.
"""
function Base.show(io::IO, ::MIME"text/plain", ms::Array{Metabolite,1})
    println(io, "Metabolite set of length: ", length(ms))
end
