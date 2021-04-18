
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
