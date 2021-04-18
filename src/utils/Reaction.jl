"""
    check_duplicate_reaction(rxn::Reaction, rxns::Dict{String, Reaction})

Check if `rxn` already exists in `rxns` but has another id.
Looks through all the reaction equations of `rxns` and compares metabolite `id`s 
and their stoichiometric coefficients to those of `rxn`.
If `rxn` has the same reaction equation as another reaction in `rxns`, the return true and id.
Otherwise return false and "".

See also: [`is_mass_balanced`](@ref)
"""
function check_duplicate_reaction(crxn::Reaction, rxns::OrderedDict{String, Reaction})
    for (k, rxn) in rxns
        if rxn.id != crxn.id # skip if same ID
            rxn_checker = true
            for (kk, vv) in rxn.metabolites # get reaction stoich
                if get(crxn.metabolites, kk, 0) != vv # if at least one stoich doesn't match
                    reaction_checker = false
                    break
                end
            end
            if reaction_checker
                return true, k
            end
        end
    end
    return false, ""
end

"""
    check_duplicate_annotations(rxn::Reaction, rxns::Dict{String, Reaction})

Determine if a `rxn` is has overlapping annotations in `rxns`.
The annotations checked are: ["bigg.reaction", "biocyc", "ec-code", "kegg.reaction", "metanetx.reaction", "rhea", "sabiork", "seed.reaction"].
Return true and the `id` of the first hit, otherwise false and "".
"""
function check_duplicate_annotations(crxn::Reaction, rxns::OrderedDict{String, Reaction})::Tuple{Bool, String}
    inspect_annotations = [
        "bigg.reaction",
        "biocyc",
        "ec-code",
        "kegg.reaction",
        "metanetx.reaction",
        "rhea",
        "sabiork",
        "seed.reaction",
    ]
    for (k, rxn) in rxns
        for anno in inspect_annotations
            if length(intersect(get(crxn.annotation, anno, ["c1"]), get(rxn.annotation, anno, ["c2"]))) != 0
                return true, k
            end
        end
    end
    return false, ""
end

"""
    is_mass_balanced(rxn::Reaction)

Checks if `rxn` is atom balanced. Returns a boolean for whether the reaction is balanced,
and the associated balance of atoms for convenience (useful if not balanced).

See also: [`get_atoms`](@ref), [`check_duplicate_reaction`](@ref)
"""
function is_mass_balanced(rxn::Reaction)
    atom_balances = Dict{String,Float64}() # float here because stoichiometry is not Int
    for (met, stoich) in rxn.metabolites
        atoms = get_atoms(met)
        isempty(atoms) && continue # ignore blanks
        for (k, v) in atoms
            atom_balances[k] = get(atom_balances, k, 0) + v * stoich
        end
    end

    return all(sum(values(atom_balances)) == 0), atom_balances
end

"""
    is_boundary(rxn::Reaction)

Return true if reaction is a boundary reaction, otherwise return false.
Checks if boundary by inspecting number of metabolites in reaction equation. 
Boundary reactions have only one metabolite, e.g. an exchange reaction, or a sink/demand reaction. 
"""
function is_boundary(rxn::Reaction)::Bool
    length(keys(rxn)) == 1 ? true : false
end
