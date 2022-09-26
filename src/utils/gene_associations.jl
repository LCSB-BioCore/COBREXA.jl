
"""
$(TYPEDSIGNATURES)

Parse `SBML.GeneProductAssociation` structure to the simpler GeneAssociation.
The input must be (implicitly) in a positive DNF.
"""
function _parse_grr(gpa::SBML.GeneProductAssociation)::GeneAssociation
    parse_ref(x) =
        typeof(x) == SBML.GPARef ? [x.gene_product] :
        begin
            @_models_log @warn "Could not parse a part of gene association, ignoring: $x"
            String[]
        end
    parse_and(x) =
        typeof(x) == SBML.GPAAnd ? vcat([parse_and(i) for i in x.terms]...) : parse_ref(x)
    parse_or(x) =
        typeof(x) == SBML.GPAOr ? vcat([parse_or(i) for i in x.terms]...) : [parse_and(x)]
    return parse_or(gpa)
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
$(TYPEDSIGNATURES)

Internal helper for parsing the string GRRs into SBML data structures. More
general than [`_parse_grr`](@ref).
"""
function _parse_grr_to_sbml(str::String)::Maybe{SBML.GeneProductAssociation}
    s = str
    toks = String[]
    m = Nothing
    while !isnothing(
        begin
            m = match(r"( +|[a-zA-Z0-9_-]+|[^ a-zA-Z0-9_()-]+|[(]|[)])(.*)", s)
        end,
    )
        tok = strip(m.captures[1])
        !isempty(tok) && push!(toks, tok)
        s = m.captures[2]
    end

    fail() = throw(DomainError(str, "Could not parse GRR"))

    # shunting yard
    ops = Symbol[]
    vals = SBML.GeneProductAssociation[]
    fold(sym, op) =
        while !isempty(ops) && last(ops) == sym
            r = pop!(vals)
            l = pop!(vals)
            pop!(ops)
            push!(vals, op([l, r]))
        end
    for tok in toks
        if tok in ["and", "AND", "&", "&&"]
            push!(ops, :and)
        elseif tok in ["or", "OR", "|", "||"]
            fold(:and, SBML.GPAAnd)
            push!(ops, :or)
        elseif tok == "("
            push!(ops, :paren)
        elseif tok == ")"
            fold(:and, SBML.GPAAnd)
            fold(:or, SBML.GPAOr)
            if isempty(ops) || last(ops) != :paren
                fail()
            else
                pop!(ops)
            end
        else
            push!(vals, SBML.GPARef(tok))
        end
    end

    fold(:and, SBML.GPAAnd)
    fold(:or, SBML.GPAOr)

    if !isempty(ops) || length(vals) > 1
        fail()
    end

    if isempty(vals)
        nothing
    else
        first(vals)
    end
end

"""
$(TYPEDSIGNATURES)

Converts a nested string gene reaction array  back into a gene reaction rule
string.

# Example
```
julia> _unparse_grr(String, [["YIL010W", "YLR043C"], ["YIL010W", "YGR209C"]])
"(YIL010W and YLR043C) or (YIL010W and YGR209C)"
```
"""
function _unparse_grr(::Type{String}, grr::GeneAssociation)::String
    grr_strings = String[]
    for gr in grr
        push!(grr_strings, "(" * join([g for g in gr], " and ") * ")")
    end
    grr_string = join(grr_strings, " or ")
    return grr_string
end
