"""
    check_duplicate_annotations(met::Metabolite, mets::OrderedDict{String, Metabolite}; inspect_annotations=_constants.metabolite_annotation_checks)

Check if a metabolite `met` has overlapping annotations with metabolites in `mets`.
The annotations checked are listed in `COBREXA._constants.metabolite_annotation_checks`.
Return id of the first hit, otherwise `nothing`.
"""
function check_duplicate_annotations(
    cmet::Metabolite,
    mets::OrderedDict{String,Metabolite};
    inspect_annotations = _constants.metabolite_annotation_checks,
)::Union{Nothing,String}
    for (k, met) in mets
        if met.compartment == cmet.compartment && k != cmet.id
            for anno in inspect_annotations
                if any(
                    in.(
                        get(met.annotations, anno, ["c1"]),
                        Ref(get(cmet.annotations, anno, ["c2"])),
                    ),
                )
                    return k
                end
            end
        end
    end
    return nothing
end

"""
    get_atoms(met::Metabolite)::MetaboliteFormula

Simple wrapper for getting the atom dictionary count out of a `Metabolite`.
"""
get_atoms(met::Metabolite)::MetaboliteFormula = _maybemap(_parse_formula, met.formula)
