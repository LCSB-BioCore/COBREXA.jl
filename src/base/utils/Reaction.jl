"""
    check_duplicate_reaction(rxn::Reaction, rxns::Dict{String, Reaction}; only_metabolites=true)

Check if `rxn` already exists in `rxns` but has another `id`.
If `only_metabolites` is `true` then only the metabolite `id`s are checked.
Otherwise, compares metabolite `id`s and the absolute value of their stoichiometric coefficients to those of `rxn`.
If `rxn` has the same reaction equation as another reaction in `rxns`, the return the `id`.
Otherwise return `nothing`.

See also: [`is_mass_balanced`](@ref)
"""
function check_duplicate_reaction(
    crxn::Reaction,
    rxns::OrderedDict{String,Reaction};
    only_metabolites = true,
)
    for (k, rxn) in rxns
        if rxn.id != crxn.id # skip if same ID
            if only_metabolites # only check if metabolites are the same
                if issetequal(keys(crxn.metabolites), keys(rxn.metabolites))
                    return k
                end
            else # also check the stoichiometric coefficients
                reaction_checker = true
                for (kk, vv) in rxn.metabolites # get reaction stoich
                    if abs(get(crxn.metabolites, kk, 0)) != abs(vv) # if at least one stoich doesn't match
                        reaction_checker = false
                        break
                    end
                end
                if reaction_checker &&
                   issetequal(keys(crxn.metabolites), keys(rxn.metabolites))
                    return k
                end
            end
        end
    end
    return nothing
end

"""
    check_duplicate_annotations(rxn::Reaction, rxns::OrderedDict{String, Reaction}; inspect_annotations=_constants.reaction_annotation_checks)

Determine if a `rxn` has overlapping annotations in `rxns`.
The annotations checked are listed in `COBREXA._constants.reaction_annotation_checks`.
Return the `id` of the first hit, otherwise `nothing`.
"""
function check_duplicate_annotations(
    crxn::Reaction,
    rxns::OrderedDict{String,Reaction};
    inspect_annotations = _constants.reaction_annotation_checks,
)::Union{Nothing,String}
    for (k, rxn) in rxns
        if k != crxn.id
            for anno in inspect_annotations
                if any(
                    occursin.(
                        get(crxn.annotations, anno, ["c1"]),
                        get(rxn.annotations, anno, ["c2"]),
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
    is_boundary(rxn::Reaction)

Return true if reaction is a boundary reaction, otherwise return false.
Checks if boundary by inspecting number of metabolites in reaction equation. 
Boundary reactions have only one metabolite, e.g. an exchange reaction, or a sink/demand reaction. 
"""
function is_boundary(rxn::Reaction)::Bool
    length(keys(rxn.metabolites)) == 1 ? true : false
end

"""
    is_mass_balanced(rxn::Reaction, model::StandardModel)

Checks if `rxn` is atom balanced. Returns a boolean for whether the reaction is balanced,
and the associated balance of atoms for convenience (useful if not balanced).

See also: [`get_atoms`](@ref), [`check_duplicate_reaction`](@ref)
"""
function is_mass_balanced(rxn::Reaction, model::StandardModel)
    atom_balances = Dict{String,Float64}() # float here because stoichiometry is not Int
    for (met, stoich) in rxn.metabolites
        atoms = get_atoms(model.metabolites[met])
        isempty(atoms) && continue # ignore blanks
        for (k, v) in atoms
            atom_balances[k] = get(atom_balances, k, 0) + v * stoich
        end
    end

    return all(sum(values(atom_balances)) == 0), atom_balances
end
