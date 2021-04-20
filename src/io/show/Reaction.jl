
"""
Pretty printing of reaction::Reaction.
"""
function Base.show(io::IO, ::MIME"text/plain", r::Reaction)
    if r.ub > 0.0 && r.lb < 0.0
        arrow = " ⟷  "
    elseif r.ub <= 0.0 && r.lb < 0.0
        arrow = " ⟵  "
    elseif r.ub > 0.0 && r.lb >= 0.0
        arrow = " ⟶  "
    else
        arrow = " →∣←  " # blocked reaction
    end
    substrates = String[]
    products = String[]
    for (k, v) in r.metabolites
        if v < 0.0
            push!(substrates, string(abs(v)) * " " * k.id)
        else
            push!(products, string(abs(v)) * " " * k.id)
        end
    end
    isempty(substrates) && (substrates = "∅")
    isempty(products) && (products = "∅")
    req_str = ""
    if length(substrates) > 5 && length(products) > 5
        sp = substrates[1] * " + ... + " * substrates[end]
        pp = products[1] * " + ... + " * products[end]
        req_str = sp * arrow * pp
    elseif length(substrates) > 5
        sp = substrates[1] * " + ... + " * substrates[end]
        req_str = sp * arrow * join(products, " + ")
    elseif length(products) > 5
        pp = products[1] * " + ... + " * products[end]
        req_str = join(substrates, " + ") * arrow * pp
    else
        req_str = join(substrates, " + ") * arrow * join(products, " + ") 
    end

    grr_strings = String[]
    for gr in r.grr
        push!(grr_strings, "(" * join([g.id for g in gr], " and ") * ")")
    end
    grr_string = join(grr_strings, " or ")
    (isnothing(grr_string) || grr_string == "") && (grr_string = "")
    
    _pretty_print(io, "Reaction ID: ", r.id)
    _pretty_print(io, "Name: ", r.name)
    _pretty_print(io, "Reaction equation: ", req_str)
    _pretty_print(io, "Lower bound: ", string(r.lb))
    _pretty_print(io, "Upper bound: ", string(r.ub))
    _pretty_print(io, "Subsystem: ", r.subsystem)
    _pretty_print(io, "Gene reaction rule: ", grr_string)
    _pretty_print(io, "Notes: ", r.notes)
    _pretty_print(io, "Annotation: ", r.annotation)
    _pretty_print(io, "Fields: ", join([string(x) for x in fieldnames(Reaction)], ", "))
end

"""
Pretty printing of reactions::Vector{Reaction}.
"""
function Base.show(io::IO, ::MIME"text/plain", rs::Vector{Reaction})
    _pretty_print(io, "Reaction vector of length: : ", string(length(rs)))
    _pretty_print(io, "Each reaction has fields: ", join([string(x) for x in fieldnames(Reaction)],", "))
end
