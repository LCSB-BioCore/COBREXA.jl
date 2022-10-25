"""
$(TYPEDSIGNATURES)

Check if `rxn` already exists in `rxns` but has another `id`.
If `only_metabolites` is `true` then only the metabolite `id`s are checked.
Otherwise, compares metabolite `id`s and the absolute value of their stoichiometric coefficients to those of `rxn`.
If `rxn` has the same reaction equation as another reaction in `rxns`, the return the `id`.
Otherwise return `nothing`.

See also: [`reaction_mass_balanced`](@ref)
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
$(TYPEDSIGNATURES)

Return true if the reaction denoted by `rxn_dict` is a boundary reaction, otherwise return false.
Checks if on boundary by inspecting the number of metabolites in `rxn_dict`.
Boundary reactions have only one metabolite, e.g. an exchange reaction, or a sink/demand reaction.
"""
function is_boundary(rxn_dict::Dict{String,Float64})::Bool
    length(keys(rxn_dict)) == 1
end

is_boundary(model::AbstractMetabolicModel, rxn_id::String) =
    is_boundary(reaction_stoichiometry(model, rxn_id))

is_boundary(rxn::Reaction) = is_boundary(rxn.metabolites)

is_boundary(model::StandardModel, rxn::Reaction) = is_boundary(rxn) # for consistency with functions below

"""
$(TYPEDSIGNATURES)

Returns a dictionary mapping the stoichiometry of atoms through a single reaction. Uses the
metabolite information in `model` to determine the mass balance. Accepts a reaction
dictionary, a reaction string id or a `Reaction` as an argument for `rxn`.

See also: [`reaction_mass_balanced`](@ref)
"""
function reaction_atom_balance(model::AbstractMetabolicModel, reaction_dict::Dict{String,Float64})
    atom_balances = Dict{String,Float64}()
    for (met, stoich_rxn) in reaction_dict
        adict = metabolite_formula(model, met)
        isnothing(adict) &&
            throw(ErrorException("Metabolite $met does not have a formula assigned to it."))
        for (atom, stoich_molecule) in adict
            atom_balances[atom] =
                get(atom_balances, atom, 0.0) + stoich_rxn * stoich_molecule
        end
    end
    return atom_balances
end

function reaction_atom_balance(model::AbstractMetabolicModel, rxn_id::String)
    reaction_atom_balance(model, reaction_stoichiometry(model, rxn_id))
end

reaction_atom_balance(model::StandardModel, rxn::Reaction) =
    reaction_atom_balance(model, rxn.id)

"""
$(TYPEDSIGNATURES)

Checks if `rxn` is atom balanced. Returns a boolean for whether the reaction is balanced,
and the associated balance of atoms for convenience (useful if not balanced). Calls
`reaction_atom_balance` internally.

See also: [`check_duplicate_reaction`](@ref), [`reaction_atom_balance`](@ref)
"""
reaction_mass_balanced(model::StandardModel, rxn_id::String) =
    all(values(reaction_atom_balance(model, rxn_id)) .== 0)

reaction_mass_balanced(model::StandardModel, rxn::Reaction) =
    reaction_mass_balanced(model, rxn.id)

reaction_mass_balanced(model::StandardModel, reaction_dict::Dict{String,Float64}) =
    all(values(reaction_atom_balance(model, reaction_dict)) .== 0)

"""
$(TYPEDSIGNATURES)

Return the reaction equation as a string. The metabolite strings can be manipulated by
setting `format_id`.

# Example
```
julia> req = Dict("coa_c" => -1, "for_c" => 1, "accoa_c" => 1, "pyr_c" => -1)
julia> stoichiometry_string(req)
"coa_c + pyr_c = for_c + accoa_c"

julia> stoichiometry_string(req; format_id = x -> x[1:end-2])
"coa + pyr = for + accoa"
```
"""
function stoichiometry_string(req; format_id = x -> x)
    count_prefix(n) = abs(n) == 1 ? "" : string(abs(n), " ")
    substrates =
        join((string(count_prefix(n), format_id(met)) for (met, n) in req if n < 0), " + ")
    products =
        join((string(count_prefix(n), format_id(met)) for (met, n) in req if n >= 0), " + ")
    return substrates * " = " * products
end

"""
$(TYPEDSIGNATURES)

Alternative of [`stoichiometry_string`](@ref) take takes a `Reaction` as an argument.
"""
stoichiometry_string(rxn::Reaction; kwargs...) =
    stoichiometry_string(rxn.metabolites; kwargs...)
