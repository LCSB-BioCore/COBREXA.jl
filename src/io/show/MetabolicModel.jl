
"""
Pretty printing of model::MetabolicModel using generic interface.
"""
function Base.show(io::IO, ::MIME"text/plain", m::MetabolicModel)
    println(Crayon(foreground=:cyan), "Metabolic model of type $(typeof(m))")
    
    print(Crayon(foreground=:blue), "Model ID: ")
    print(Crayon(foreground=:magenta), "Use generic function to get ID")

    if prod(size(stoichiometry(m))) < 5_000_000
        println(Crayon(foreground=:blue), stoichiometry(m))
    else # too big to display nicely
        println(Crayon(foreground=:blue), "S × v = b")    
    end
    
    print(Crayon(foreground=:blue), "Number of reactions: ")
    println(Crayon(foreground=:magenta), n_reactions(m))
    
    print(Crayon(foreground=:blue), "Number of metabolites: ")
    println(Crayon(foreground=:magenta),  n_metabolites(m))

    print(Crayon(foreground=:blue), "Field names associated with this model include: ")
    print(Crayon(foreground=:magenta),join([string(x) for x in fieldnames(typeof(m))], ", "))
end
