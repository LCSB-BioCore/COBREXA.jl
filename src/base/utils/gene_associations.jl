
"""
    _parse_grr(gpa::SBML.GeneProductAssociation)::GeneAssociation

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
    _unparse_grr(
        ::Type{SBML.GeneProductAssociation},
        x::GeneAssociation,
    )::SBML.GeneAssociation

Convert a GeneAssociation to the corresponding `SBML.jl` structure.
"""
function _unparse_grr(
    ::Type{SBML.GeneProductAssociation},
    x::GeneAssociation,
)::SBML.GeneAssociation
    SBML.GPAOr([SBML.GPAAnd([SBML.GPARef(j) for j in i]) for i in x])
end

"""
    _parse_grr(s::String)::GeneAssociation

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
    _parse_grr_to_sbml(str::String)::Maybe{SBML.GeneProductAssociation}

Internal helper for parsing the string GRRs into SBML data structures. More
general than [`_parse_grr`](@ref).
"""
function _parse_grr_to_sbml(str::String)::Maybe{SBML.GeneProductAssociation}
    s = str
    toks = String[]
    m = Nothing
    while !isnothing(begin
        m = match(r"( +|[a-zA-Z0-9]+|[^ a-zA-Z0-9()]+|[(]|[)])(.*)", s)
    end)
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
    unparse_grr(grr::Vector{Vector{Gene}}

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

"""
    enzyme_availability(
        grr::GeneAssociation,
        gene_availability::Function;
        or_function,
        and_function,
    )

Compute the amount of enzyme available for catalyzing a reaction from
[`GeneAssociation`](@ref) `grr` and amount of genes available. The function
will use `gene_available` to query for gene identifiers strings in order to
receive the gene availability values, combine these using `and_function` and
`or_function` for the levels of DNF grr, and return the result. Other overloads
of the method provide simplified interface for querying dictionaries, and
pre-set functions for combining the expressions.
"""
enzyme_availability(
    grr::GeneAssociation,
    gene_availability::Function;
    or_function,
    and_function,
)::Float64 = or_function(map(gs -> and_function(map(gene_availability, gs)), grr))

"""
    enzyme_availability(
        grr::GeneAssociation,
        gene_availability::Dict{String,A};
        kwargs...,
    )

Overload of [`enzyme_availability`](@ref) that takes the amounts of available
genes from a dictionary. Error is thrown in case of genes missing in the
dictionary.
"""
enzyme_availability(
    grr::GeneAssociation,
    gene_availability::Dict{String,A};
    kwargs...,
) where {A} = enzyme_availability(grr, g -> gene_availability[g]; kwargs...)

"""
    enzyme_availability(
        grr::GeneAssociation,
        gene_availability::Dict{String,A},
        default::A;
        kwargs...,
    )

Overload of [`enzyme_availability`](@ref) that takes the amounts of available
genes from a dictionary. The `default` is used for values missing in the
dictionary.
"""
enzyme_availability(
    grr::GeneAssociation,
    gene_availability::Dict{String,A},
    default::A;
    kwargs...,
) where {A} = enzyme_availability(grr, g -> get(gene_availability, g, default); kwargs...)

"""
    enzyme_availability_probabilistic(args...; kwargs...)

Variant of [`enzyme_availability`](@ref) that sets up `and_function` and
`or_function` so that the amounts of genes are interpreted as probabilities.
Forwards all arguments to [`enzyme_availability`](@ref).
"""
enzyme_availability_probabilistic(args...; kwargs...) = enzyme_availability(
    args...;
    or_function = vs -> 1 - prod(1 .- vs),
    and_function = prod,
    kwargs...,
)

"""
    enzyme_availability_eflux(args...; kwargs...)

Variant of [`enzyme_availability`](@ref) that sets up `and_function` and
`or_function` so that the relative amounts of genes material are interpreted
similarly as in E-Flux algorithm, taking a minimum of each gene group, then
adding the results.  Forwards all arguments to [`enzyme_availability`](@ref).

For details, see: Colijn, Caroline, et al. "Interpreting expression data with
metabolic flux models: predicting Mycobacterium tuberculosis mycolic acid
production." PLoS computational biology 5.8 (2009): e1000489.
"""
enzyme_availability_eflux(args...; kwargs...) =
    enzyme_availability(args...; or_function = sum, and_function = min, kwargs...)
