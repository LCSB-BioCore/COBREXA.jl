"""
    getindex(rxns::Vector{Reaction}, rxn::Reaction)

Get the index of a reaction `rxn` in an array of reactions `rxns`, based in `id`.
Return -1 if no matches found.

Typically used, `index = rxns[rxn]`.
"""
function Base.getindex(rxns::Vector{Reaction}, rxn::Reaction)
    for i in eachindex(rxns)
        if rxns[i].id == rxn.id
            return i
        end
    end
    return -1
end

"""
    findfirst(rxns::Vector{Reaction}, rxnid::String)

Return the reaction with `rxnid` in `rxns` or else `nothing`.

Typically used: `findfirst(model.rxns, rxnid)`.
"""
function Base.findfirst(rxns::Vector{Reaction}, rxnid::String)
    for i in eachindex(rxns)
        if rxns[i].id == rxnid
            return rxns[i]
        end
    end
    return nothing
end

"""
    check_duplicate_reaction(rxns::Vector{Reaction}, rxn::Reaction)

Check if `rxn` already exists in `rxns` but has another id.
Looks through all the reaction equations of `rxns` and compares metabolite `id`s and their stoichiometric coefficients to those of `rxn`.
If `rxn` has the same reaction equation as another reaction in `rxns`, the return true and the index of the first match.

See also: [`is_mass_balanced`](@ref)
"""
function check_duplicate_reaction(rxns::Vector{Reaction}, crxn::Reaction)
    ceq = Dict{String,Float64}(k.id => v for (k, v) in crxn.metabolites)
    for rxn in rxns
        if rxn.id != crxn.id
            req = Dict{String,Float64}(k.id => v for (k, v) in rxn.metabolites)
            if isempty(setdiff(collect(keys(ceq)), collect(keys(req)))) &&
               isempty(setdiff(collect(keys(req)), collect(keys(ceq)))) # same metabolites
                all([req[k] == v for (k, v) in ceq]) && (return true, rxns[rxn])
            end
        else
            return true, rxns[rxn]
        end
    end
    return false, -1
end

"""
    check_duplicate_annotations(rxns::Vector{Reaction}, rxn::Gene)

Determine if a `rxn` is has overlapping annotations in `rxns`.
The annotations checked are: ["bigg.reaction", "biocyc", "ec-code", "kegg.reaction", "metanetx.reaction", "rhea", "sabiork", "seed.reaction"].
Return true and the index of the first hit, otherwise false and -1.
"""
function check_duplicate_annotations(rxns::Vector{Reaction}, crxn::Reaction)
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
    for rxn in rxns
        for anno in inspect_annotations
            if any([
                x in get(crxn.annotation, anno, ["c1"]) for
                x in get(rxn.annotation, anno, ["c2"])
            ])
                return true, rxns[rxn]
            end
        end
    end
    return false, -1
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
