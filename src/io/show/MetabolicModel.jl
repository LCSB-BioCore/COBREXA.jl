
"""
Pretty printing of model::MetabolicModel using generic interface.
"""
function Base.show(io::IO, ::MIME"text/plain", m::MetabolicModel)
    println(Crayon(foreground=:cyan), "Metabolic model: [INSERT ID]")
    println(Crayon(foreground=:blue), stoichiometry(m))
    
    print(Crayon(foreground=:blue), "Number of reactions: ")
    println(Crayon(foreground=:magenta), n_reactions(m))
    
    print(Crayon(foreground=:blue), "Number of metabolites: ")
    println(Crayon(foreground=:magenta),  n_reactions(m))

    print(Crayon(foreground=:blue), "Field names associated with this model include: ")
    print(Crayon(foreground=:magenta),join([string(x) for x in fieldnames(typeof(m))], ", "))
end
