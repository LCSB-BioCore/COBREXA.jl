
"""
Pretty printing of model::StandardModel.
"""
function Base.show(io::IO, ::MIME"text/plain", m::StandardModel)
    println(
        io,
        "Constraint based model: ",
        m.id,
        "\n",
        "Number of reactions: ",
        length(m.reactions),
        "\n",
        "Number of metabolites: ",
        length(m.metabolites),
        "\n",
        "Number of genes: ",
        length(m.genes),
    )
end
