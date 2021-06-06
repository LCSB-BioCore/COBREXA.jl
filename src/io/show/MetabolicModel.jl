"""
Pretty printing of everything metabolic-modelish.
"""
function Base.show(io::IO, ::MIME"text/plain", m::MetabolicModel)
    _pretty_print_keyvals(io, "", "Metabolic model of type $(typeof(m))")
    if prod(size(stoichiometry(m))) < _constants.default_stoichiometrix_view_size
        println(io, stoichiometry(m))
    else # too big to display nicely
        println(io, "S = [...]")
    end
    _pretty_print_keyvals(io, "Number of reactions: ", string(n_reactions(m)))
    _pretty_print_keyvals(io, "Number of metabolites: ", string(n_metabolites(m)))
end
