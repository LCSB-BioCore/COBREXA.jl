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
            push!(substrates, string(abs(v)) * " " * k)
        else
            push!(products, string(abs(v)) * " " * k)
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
        push!(grr_strings, "(" * join([g for g in gr], " and ") * ")")
    end
    grr_string = join(grr_strings, " or ")
    (isnothing(grr_string) || grr_string == "") && (grr_string = "")

    for fname in fieldnames(Reaction)
        if fname == :metabolites
            _print_color(io, "Reaction.$(string(fname)): ", req_str)
        elseif fname == :grr
            _print_color(io, "Reaction.$(string(fname)): ", grr_string)
        elseif fname in (:lb, :ub, :objective_coefficient)
            _print_color(io, "Reaction.$(string(fname)): ", string(getfield(r, fname)))
        else
            _print_color(io, "Reaction.$(string(fname)): ", getfield(r, fname))
        end
    end
end
