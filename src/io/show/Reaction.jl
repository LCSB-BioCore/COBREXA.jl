
"""
    _pretty_substances(ss::Vector{String})::String

Nicely format a substance list.
"""
function _pretty_substances(ss::Vector{String})::String
    if isempty(ss)
        "∅"
    elseif length(ss) > 5
        join([ss[1], ss[2], "...", ss[end-1], ss[end]], " + ")
    else
        join(ss, " + ")
    end
end

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
    substrates =
        ["$(-v) $k" for (k, v) in Iterators.filter(((_, v)::Pair -> v < 0), r.metabolites)]
    products =
        ["$v $k" for (k, v) in Iterators.filter(((_, v)::Pair -> v >= 0), r.metabolites)]

    for fname in fieldnames(Reaction)
        if fname == :metabolites
            _print_with_colors(
                io,
                "Reaction.$(string(fname)): ",
                _pretty_substances(substrates) * arrow * _pretty_substances(products),
            )
        elseif fname == :grr
            _print_with_colors(
                io,
                "Reaction.$(string(fname)): ",
                _maybemap(x -> _unparse_grr(String, x), r.grr),
            )
        elseif fname in (:lb, :ub, :objective_coefficient)
            _print_with_colors(
                io,
                "Reaction.$(string(fname)): ",
                string(getfield(r, fname)),
            )
        else
            _print_with_colors(io, "Reaction.$(string(fname)): ", getfield(r, fname))
        end
    end
end
