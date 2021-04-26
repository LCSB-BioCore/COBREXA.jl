
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

    _print_color(io, "Reaction ID: ", r.id)
    _print_color(io, "Name: ", r.name)
    _print_color(io, "Reaction equation: ", req_str)
    _print_color(io, "Lower bound: ", string(r.lb))
    _print_color(io, "Upper bound: ", string(r.ub))
    _print_color(io, "Subsystem: ", r.subsystem)
    _print_color(io, "Gene reaction rule: ", grr_string)
    _print_color(io, "Notes: ", r.notes)
    _print_color(io, "Annotation: ", r.annotation)
    _print_color(io, "Fields: ", join([string(x) for x in fieldnames(Reaction)], ", "))
end

"""
Pretty printing of reactions::Vector{Reaction}.
"""
function Base.show(io::IO, ::MIME"text/plain", rs::Vector{Reaction})
    _print_color(io, "Reaction vector of length: : ", string(length(rs)))
    _print_color(
        io,
        "Each reaction has fields: ",
        join([string(x) for x in fieldnames(Reaction)], ", "),
    )
end
