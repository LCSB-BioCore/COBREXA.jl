
"""
Pretty printing of model::MetabolicModel using generic interface.
"""
function Base.show(io::IO, ::MIME"text/plain", m::MetabolicModel)
    println(io, Crayon(foreground=:cyan), "Metabolic model of type $(typeof(m))")
    
    print(io, Crayon(foreground=:blue), "Model ID: ")
    print(io, Crayon(foreground=:magenta), "Use generic function to get ID")

    if prod(size(stoichiometry(m))) < 5_000_000
        println(io, Crayon(foreground=:blue), stoichiometry(m))
    else # too big to display nicely
        println(io, Crayon(foreground=:blue), "S × v = b")    
    end
    
    print(io, Crayon(foreground=:blue), "Number of reactions: ")
    println(io, Crayon(foreground=:magenta), n_reactions(m))
    
    print(io, Crayon(foreground=:blue), "Number of metabolites: ")
    println(io, Crayon(foreground=:magenta),  n_metabolites(m))

    print(io, Crayon(foreground=:blue), "Field names associated with this model include: ")
    print(io, Crayon(foreground=:magenta),join([string(x) for x in fieldnames(typeof(m))], ", "))
end
