
"""
Pretty printing of model::MetabolicModel using generic interface.
"""
function Base.show(io::IO, ::MIME"text/plain", m::MetabolicModel)
    println(io, Crayon(foreground=:cyan), "Metabolic model of type $(typeof(m))")
    _pretty_print(io, "Model ID: ", "Use generic function to get ID")

    if prod(size(stoichiometry(m))) < 5_000_000
        println(io, Crayon(foreground=:blue), stoichiometry(m))
    else # too big to display nicely
        println(io, Crayon(foreground=:blue), "S × v = b")    
    end
    
    _pretty_print(io, "Number of reactions: ", string(n_reactions(m)))
    _pretty_print(io, "Number of metabolites: ", string(n_metabolites(m)))
    _pretty_print(io, "Field names associated with this model include: ", join([string(x) for x in fieldnames(typeof(m))], ", "))
end
