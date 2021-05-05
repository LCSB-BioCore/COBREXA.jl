
"""
    _parse_formula(f::String)::MetaboliteFormula

Parse a formula in format `C2H6O` into a [`MetaboliteFormula`](@ref), which is
basically a dictionary of atom counts in the molecule.
"""
function _parse_formula(f::String)::MetaboliteFormula
    res = Dict{String,Int}()
    pattern = @r_str "([A-Z][a-z]*)([1-9][0-9]*)?"

    for m in eachmatch(pattern, f)
        res[m.captures[1]] = isnothing(m.captures[2]) ? 1 : parse(Int, m.captures[2])
    end

    return res
end

"""
    _unparse_formula(f::MetaboliteFormula)::String

Format [`MetaboliteFormula`](@ref) to `String`.
"""
function _unparse_formula(f::MetaboliteFormula)::String
    return join(["$elem$n" for (elem, n) in f])
end
