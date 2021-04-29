
"""
    _formula_to_atoms(f::String)::Dict{String,Int}
"""
function _formula_to_atoms(f::String)::MetaboliteFormula
    res = Dict{String,Int}()
    pattern = @r_str "([A-Z][a-z]*)([1-9][0-9]*)?"

    for m in eachmatch(pattern, f)
        res[m.captures[1]] = isnothing(m.captures[2]) ? 1 : parse(Int, m.captures[2])
    end

    return res
end

"""
    _atoms_to_formula(f::String)::Dict{String,Int}
"""
function _atoms_to_formula(f::MetaboliteFormula)::String
    return join(["$elem$n" for (elem, n) in f])
end

function _atoms_to_formula(_::Nothing)::Nothing
    return nothing
end
