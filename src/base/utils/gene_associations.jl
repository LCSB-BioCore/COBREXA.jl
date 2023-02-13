
"""
$(TYPEDSIGNATURES)

A helper for producing predictable unique sequences. Might be faster if
compacting would be done directly in sort().
"""
function _sortunique(x)
    o = collect(x)
    sort!(o)
    put = prevind(o, firstindex(o))
    for i in eachindex(o)
        if put >= firstindex(o) && o[i] == o[put]
            # we already have this one
            continue
        else
            put = nextind(o, put)
            if put != i
                o[put] = o[i]
            end
        end
    end
    o[begin:put]
end

"""
$(TYPEDSIGNATURES)

Parse `SBML.GeneProductAssociation` structure and convert it to a strictly
positive DNF [`GeneAssociation`](@ref). Negation (`SBML.GPANot`) is not
supported.
"""
function _parse_grr(gpa::SBML.GeneProductAssociation)::GeneAssociation

    function fold_and(dnfs::Vector{Vector{Vector{String}}})::Vector{Vector{String}}
        if isempty(dnfs)
            [String[]]
        else
            _sortunique(
                _sortunique(String[l; r]) for l in dnfs[1] for r in fold_and(dnfs[2:end])
            )
        end
    end

    dnf(x::SBML.GPARef) = [[x.gene_product]]
    dnf(x::SBML.GPAOr) = _sortunique(vcat(dnf.(x.terms)...))
    dnf(x::SBML.GPAAnd) = fold_and(dnf.(x.terms))
    dnf(x) = throw(
        DomainError(
            x,
            "unsupported gene product association contents of type $(typeof(x))",
        ),
    )
    return dnf(gpa)
end

"""
$(TYPEDSIGNATURES)

Convert a GeneAssociation to the corresponding `SBML.jl` structure.
"""
function _unparse_grr(
    ::Type{SBML.GeneProductAssociation},
    x::GeneAssociation,
)::SBML.GeneProductAssociation
    SBML.GPAOr([SBML.GPAAnd([SBML.GPARef(j) for j in i]) for i in x])
end

"""
$(TYPEDSIGNATURES)

Parse a DNF gene association rule in format `(YIL010W and YLR043C) or (YIL010W
and YGR209C)` to `GeneAssociation. Also accepts `OR`, `|`, `||`, `AND`, `&`,
and `&&`.

# Example
```
julia> _parse_grr("(YIL010W and YLR043C) or (YIL010W and YGR209C)")
2-element Array{Array{String,1},1}:
 ["YIL010W", "YLR043C"]
 ["YIL010W", "YGR209C"]
```
"""
_parse_grr(s::String)::Maybe{GeneAssociation} = _maybemap(_parse_grr, _parse_grr_to_sbml(s))

"""
PikaParser grammar for stringy GRR expressions.
"""
const _grr_grammar = begin
    # characters that typically form the identifiers
    isident(x::Char) =
        isletter(x) ||
        isdigit(x) ||
        x == '_' ||
        x == '-' ||
        x == ':' ||
        x == '.' ||
        x == '\'' ||
        x == '[' ||
        x == ']' ||
        x == '\x03' # a very ugly exception for badly parsed MAT files

    # scanner helpers
    eat(p) = m -> begin
        last = 0
        for i in eachindex(m)
            p(m[i]) || break
            last = i
        end
        last
    end

    # eat one of keywords
    kws(w...) = m -> begin
        last = eat(isident)(m)
        m[begin:last] in w ? last : 0
    end

    PP.make_grammar(
        [:expr],
        PP.flatten(
            Dict(
                :space => PP.first(PP.scan(eat(isspace)), PP.epsilon),
                :id => PP.scan(eat(isident)),
                :orop =>
                    PP.first(PP.tokens("||"), PP.token('|'), PP.scan(kws("OR", "or"))),
                :andop => PP.first(
                    PP.tokens("&&"),
                    PP.token('&'),
                    PP.scan(kws("AND", "and")),
                ),
                :expr => PP.seq(:space, :orexpr, :space, PP.end_of_input),
                :orexpr => PP.first(
                    :or => PP.seq(:andexpr, :space, :orop, :space, :orexpr),
                    :andexpr,
                ),
                :andexpr => PP.first(
                    :and => PP.seq(:baseexpr, :space, :andop, :space, :andexpr),
                    :baseexpr,
                ),
                :baseexpr => PP.first(
                    :id,
                    :parenexpr => PP.seq(
                        PP.token('('),
                        :space,
                        :orexpr,
                        :space,
                        PP.token(')'),
                    ),
                ),
            ),
            Char,
        ),
    )
end

_grr_grammar_open(m, _) =
    m.rule == :expr ? Bool[0, 1, 0, 0] :
    m.rule == :parenexpr ? Bool[0, 0, 1, 0, 0] :
    m.rule in [:or, :and] ? Bool[1, 0, 0, 0, 1] :
    m.rule in [:andexpr, :orexpr, :notexpr, :baseexpr] ? Bool[1] :
    (false for _ in m.submatches)

_grr_grammar_fold(m, _, subvals) =
    m.rule == :id ? SBML.GPARef(m.view) :
    m.rule == :and ? SBML.GPAAnd([subvals[1], subvals[5]]) :
    m.rule == :or ? SBML.GPAOr([subvals[1], subvals[5]]) :
    m.rule == :parenexpr ? subvals[3] :
    m.rule == :expr ? subvals[2] : isempty(subvals) ? nothing : subvals[1]

"""
$(TYPEDSIGNATURES)

Internal helper for parsing the string GRRs into SBML data structures. More
general than [`_parse_grr`](@ref).
"""
function _parse_grr_to_sbml(str::String)::Maybe{SBML.GeneProductAssociation}
    all(isspace, str) && return nothing
    tree = PP.parse_lex(_grr_grammar, str)
    match = PP.find_match_at!(tree, :expr, 1)
    if match > 0
        return PP.traverse_match(
            tree,
            match,
            open = _grr_grammar_open,
            fold = _grr_grammar_fold,
        )
    else
        throw(DomainError(str, "cannot parse GRR"))
    end
end

"""
$(TYPEDSIGNATURES)

Converts a nested string gene reaction array  back into a gene reaction rule
string.

# Example
```
julia> _unparse_grr(String, [["YIL010W", "YLR043C"], ["YIL010W", "YGR209C"]])
"(YIL010W && YLR043C) || (YIL010W && YGR209C)"
```
"""
function _unparse_grr(
    ::Type{String},
    grr::GeneAssociation;
    and = " && ",
    or = " || ",
)::String
    return join(("(" * join(gr, and) * ")" for gr in grr), or)
end
