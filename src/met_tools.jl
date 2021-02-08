function getatoms(met::Metabolite)
    atoms = Dict{String, Int64}()
    N = length(met.formula)

    caps = findall(x-> isuppercase(x[1]), split(met.formula, ""))
    if length(caps) == 1
        atom, count = formulapart(met.formula)
        atoms[atom] = count
    else
        pend = [caps[2:end].-1;length(met.formula)]
        for (i, j) in zip(caps, pend) 
            atom, count = formulapart(met.formula[i:j])
            atoms[atom] = count            
        end
    end
    return atoms
end

function formulapart(formula)
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