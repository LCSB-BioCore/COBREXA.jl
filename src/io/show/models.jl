"""
Pretty printing of everything metabolic-modelish.
"""
function Base.show(
    io::IO,
    ::MIME"text/plain",
    m::MetabolicModel,
)
    _print_with_colors(io, "", "Metabolic model of type $(typeof(m))")
    if prod(size(stoichiometry(m))) < 5_000_000
        println(io, Crayon(foreground = _constants.colors.payload), stoichiometry(m))
    else # too big to display nicely
        println(io, Crayon(foreground = _constants.colors.payload), "S = [...]")
    end
    _print_with_colors(io, "Number of reactions: ", string(n_reactions(m)))
    _print_with_colors(io, "Number of metabolites: ", string(n_metabolites(m)))
end
