
function _pretty_substances(ss::Vector{String})::String
    if isempty(ss)
        "∅"
    elseif length(ss) > 5
        join([ss[1], ss[2], "…", ss[end-1], ss[end]], " + ")
    else
        join(ss, " + ")
    end
end

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
    substrates = [
        "$(-v) $(k.id)" for
        (k, v) in Iterators.filter(((_, v)::Pair -> v < 0), r.metabolites)
    ]
    products = [
        "$v $(k.id)" for
        (k, v) in Iterators.filter(((_, v)::Pair -> v >= 0), r.metabolites)
    ]

    grr_strings = String[]
    for gr in r.grr
        push!(grr_strings, "(" * join([g for g in gr], " and ") * ")")
    end
    grr_string = join(grr_strings, " or ")
    (isnothing(grr_string) || grr_string == "") && (grr_string = "")

    for fname in fieldnames(Reaction)
        if fname == :metabolites
            _print_color(
                io,
                "Reaction.$(string(fname)): ",
                _format_substances(substrates) * arrow * _format_substances(products),
            )
        elseif fname == :grr
            _print_color(io, "Reaction.$(string(fname)): ", grr_string)
        elseif fname in (:lb, :ub, :objective_coefficient)
            _print_color(io, "Reaction.$(string(fname)): ", string(getfield(r, fname)))
        else
            _print_color(io, "Reaction.$(string(fname)): ", getfield(r, fname))
        end
    end
end
