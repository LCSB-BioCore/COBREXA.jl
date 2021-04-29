"""
    check_duplicate_reaction(rxn::Reaction, rxns::Dict{String, Reaction})

Check if `rxn` already exists in `rxns` but has another `id`.
Looks through all the reaction equations of `rxns` and compares metabolite `id`s 
and their stoichiometric coefficients to those of `rxn`.
If `rxn` has the same reaction equation as another reaction in `rxns`, the return the `id`.
Otherwise return `nothing`.

See also: [`is_mass_balanced`](@ref)
"""
function check_duplicate_reaction(crxn::Reaction, rxns::OrderedDict{String,Reaction})
    for (k, rxn) in rxns
        if rxn.id != crxn.id # skip if same ID
            reaction_checker = true
            for (kk, vv) in rxn.metabolites # get reaction stoich
                if get(crxn.metabolites, kk, 0) != vv # if at least one stoich doesn't match
                    reaction_checker = false
                    break
                end
            end
            if reaction_checker
                return k
            end
        end
    end
    return nothing
end

"""
    check_duplicate_annotations(rxn::Reaction, rxns::Dict{String, Reaction})

Determine if a `rxn` is has overlapping annotations in `rxns`.
The annotations checked are: ["bigg.reaction", "biocyc", "ec-code", "kegg.reaction", "metanetx.reaction", "rhea", "sabiork", "seed.reaction"].
Return true and the `id` of the first hit, otherwise false and "".
"""
function check_duplicate_annotations(
    crxn::Reaction,
    rxns::OrderedDict{String,Reaction};
    inspect_annotations = [
        "bigg.reaction",
        "biocyc",
        "ec-code",
        "kegg.reaction",
        "metanetx.reaction",
        "rhea",
        "sabiork",
        "seed.reaction",
    ],
)::Union{Nothing,String}
    for (k, rxn) in rxns
        for anno in inspect_annotations
            if length(
                intersect(
                    get(crxn.annotations, anno, ["c1"]),
                    get(rxn.annotations, anno, ["c2"]),
                ),
            ) != 0
                return k
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
