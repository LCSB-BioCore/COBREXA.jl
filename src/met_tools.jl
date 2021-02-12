"""
get_atoms(met::Metabolite)

Return a dictionary mapping the elements in a metabolite to their stoichiometric coefficients.
"""
function get_atoms(met::Metabolite)
    atoms = Dict{String, Int64}()
    N = length(met.formula)
    
    N == 0 && return atoms

    caps = findall(x-> isuppercase(x[1]), split(met.formula, ""))
    if length(caps) == 1
        atom, count = split_formula(met.formula)
        atoms[atom] = count
    else
        pend = [caps[2:end].-1;length(met.formula)]
        for (i, j) in zip(caps, pend) 
            atom, count = split_formula(met.formula[i:j])
            atoms[atom] = count
        end
    end
    return atoms
end

"""
split_formula(formula)

Split the Atom from the stoichiometric coefficient. 
E.g. C12 => C, 12 or Ca3 => Ca, 3
"""
function split_formula(formula)
    N = length(formula)
    if N > 1
        if islowercase(formula[2][1])
            if N > 2
                atom = string(formula[1:2])
                count = parse(Int64, formula[3:end])
            else
                atom = string(formula[1:2])
                count = 1
            end
        else
            atom = string(formula[1])
            count = parse(Int64, formula[2:end]) 
        end
    else
        atom = string(formula[1])
        count = 1
    end
    return atom, count
end