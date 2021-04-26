
"""
Pretty printing of model::StandardModel.
"""
function Base.show(io::IO, ::MIME"text/plain", m::StandardModel)
    println(io, "Constraint based model: $(m.id)\nNumber of reactions: $(length(m.reactions))\nNumber of metabolites: $(length(m.metabolites))\nNumber of genes: $(length(m.genes))")
end
