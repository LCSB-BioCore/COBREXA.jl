
"""
    _parse_grr(gpa::SBML.GeneProductAssociation)::GeneAssociation

Parse `SBML.GeneProductAssociation` structure to the simpler GeneAssociation.
The input must be (implicitly) in a positive DNF.
"""
function _parse_grr(gpa::SBML.GeneProductAssociation)::GeneAssociation
    parse_ref(x) = typeof(x) == SBML.GPARef ? x.gene_product : nothing
    parse_ands(x) =
        typeof(x) == SBML.GPAAnd ? [parse_ref(i) for i in x.terms] : parse_ref(x)
    parse_or(x) = typeof(x) == SBML.GPAOr ? [parse_and(i) for i in x.terms] : parse_and(x)

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
function _parse_grr(s::String)::GeneAssociation
    # first get the gene id list in string format
    gene_reaction_rules = Vector{Vector{String}}()
    or_genes = split(s, r"\s?(or|OR|(\|\|)|\|)\s?") # separate or terms
    for or_gene in or_genes
        and_genes = split(replace(or_gene, r"\(|\)" => ""), r"\s?(and|AND|(\&\&)|\&)\s?")
        push!(gene_reaction_rules, and_genes)
    end
    return gene_reaction_rules
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
