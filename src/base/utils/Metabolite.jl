"""
    get_atoms(met::Metabolite)::MetaboliteFormula

Simple wrapper for getting the atom dictionary count out of a `Metabolite`.
"""
get_atoms(met::Metabolite)::MetaboliteFormula = _maybemap(_parse_formula, met.formula)
