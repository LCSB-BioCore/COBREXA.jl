"""
    get_atoms(met::Metabolite)::MetaboliteFormula

Simple wrapper for getting the atom dictionary count out of a `Metabolite`.

See also: [`metabolite_formula`](@ref)
"""
get_atoms(met::Metabolite)::MetaboliteFormula = _maybemap(_parse_formula, met.formula)
