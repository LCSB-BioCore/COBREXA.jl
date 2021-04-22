
"""
Pretty printing of reaction::Reaction.
"""
function Base.show(io::IO, ::MIME"text/plain", r::Reaction)
    if r.ub > 0.0 && r.lb < 0.0
        arrow = " ⟷  "
    elseif r.ub == 0.0 && r.lb < 0.0
        arrow = " ⟵  "
    elseif r.ub > 0.0 && r.lb == 0.0
        arrow = " ⟶  "
    else
        arrow = " →∣←  " # blocked reaction
    end
    substrates = String[]
    products = String[]
    for (k, v) in r.metabolites
        if v < 0.0
            push!(substrates, string(abs(v)) * " " * k)
        else
            push!(products, string(abs(v)) * " " * k)
        end
    end
    isempty(substrates) && (substrates = "∅")
    isempty(products) && (products = "∅")

    println(io, "Reaction ID: ", r.id)
    println(io, "Reaction name: ", r.name)
    println(io, "Reaction subsystem: ", r.subsystem)
    if length(substrates) > 5 && length(products) > 5
        sp = substrates[1] * " + ... + " * substrates[end]
        pp = products[1] * " + ... + " * products[end]
        println(io, sp * arrow * pp)
    elseif length(substrates) > 5
        sp = substrates[1] * " + ... + " * substrates[end]
        println(io, sp * arrow * join(products, " + "))
    elseif length(products) > 5
        pp = products[1] * " + ... + " * products[end]
        println(io, join(substrates, " + ") * arrow * pp)
    else
        println(io, join(substrates, " + ") * arrow * join(products, " + "))
    end
    println(io, "Lower bound: ", r.lb)
    println(io, "Upper bound: ", r.ub)

    grr_strings = String[]
    for gr in r.grr
        push!(grr_strings, "(" * join([g for g in gr], " and ") * ")")
    end
    grr_string = join(grr_strings, " or ")
    (isnothing(grr_string) || grr_string == "") && (grr_string = "")
    println(io, "Genes: ", grr_string)
    println(io, "E.C. number: ", join(get(r.annotation, "ec-code", [""]), " or "))
end

"""
Pretty printing of reactions::Vector{Reaction}.
"""
function Base.show(io::IO, ::MIME"text/plain", rs::Vector{Reaction})
    println(io, "Reaction set of length: ", length(rs))
end
