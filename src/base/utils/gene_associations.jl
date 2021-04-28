
#TODO: convert the string-GeneAssociation formatting to this too

function _parse_grr(gpa::SBML.GeneProductAssociation)::GeneAssociation
    parse_ref(x) = typeof(x) == SBML.GPARef ? x.gene_product : nothing
    parse_ands(x) =
        typeof(x) == SBML.GPAAnd ? [parse_ref(i) for i in x.terms] : parse_ref(x)
    parse_or(x) = typeof(x) == SBML.GPAOr ? [parse_and(i) for i in x.terms] : parse_and(x)

    return parse_or(gpa)
end

function _unparse_grr(
    ::Type{SBML.GeneProductAssociation},
    x::GeneAssociation,
)::SBML.GeneAssociation
    SBML.GPAOr([SBML.GPAAnd([SBML.GPARef(j) for j in i]) for i in x])
end

"""
    parsegrr(string_rule)

Parse a gene reaction rule string `string_rule` into a nested string gene reaction rule array `Vector{Vector{String}}`.

Format: (YIL010W and YLR043C) or (YIL010W and YGR209C) where `or` can also be `OR, |, ||` and where `and` can also be `AND, &, &&`.
So `"(YIL010W and YLR043C) or (YIL010W and YGR209C)"`` becomes `[[YIL010W, YLR043C], [YIL010W, YGR209C]]`` as a nested string array.
"""
function _parse_grr(s::Maybe{String})::Maybe{GeneAssociation}
    if s == "" || isnothing(s)
        return nothing
    end
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

Converts a nested string gene reaction array  back into a gene reaction rule string.
So `[[YIL010W, YLR043C], [YIL010W, YGR209C]]`` becomes `"(YIL010W and YLR043C) or (YIL010W and YGR209C)"``.
"""
function _unparse_grr(grr::Maybe{GeneAssociation})::Maybe{String}
    isnothing(grr) && return nothing
    grr_strings = String[]
    for gr in grr
        push!(grr_strings, "(" * join([g for g in gr], " and ") * ")")
    end
    grr_string = join(grr_strings, " or ")
    return grr_string
end
