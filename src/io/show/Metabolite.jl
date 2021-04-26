
"""
Pretty printing of metabolite::Metabolite.
"""
function Base.show(io::IO, ::MIME"text/plain", m::Metabolite)
    println(io, "Metabolite ID: $(m.id)\nMetabolite name: $(m.name)\nFormula: $(m.formula)\nCharge: $(m.charge)")
end

"""
Pretty printing of metabolites::Vector{Metabolite}.
"""
function Base.show(io::IO, ::MIME"text/plain", ms::Vector{Metabolite})
    println(io, "Metabolite set of length: ", length(ms))
end
