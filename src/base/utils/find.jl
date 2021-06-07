"""
    function find_exchange_reactions(
        model::MetabolicModel;
        exclude_biomass = false,
        biomass_strings = _constants.biomass_strings,
        ex_prefixes = _constants.exchange_prefixes,
    )::Vector{String}

Return a vector of reaction ids of exchange reactions. Identifies these reactions based
on their prefixes, see `_constants.exchange_prefixes` for the default list (includes 
prefixes like "EX_", etc.). If `exclude_biomass` is true then the biomass reaction is
not returned as in this list, otherwise looks for the biomass reaction by checking if 
a reaction contains a string in the list `biomass_strings`, see `_constants.biomass_strings`
for the default list (includes strings like "biomass", etc.).
"""
function find_exchange_reactions(
    model::M;
    exclude_biomass = false,
    biomass_strings = _constants.biomass_strings,
    ex_prefixes = _constants.exchange_prefixes,
)::Vector{String} where M <: MetabolicModel
    ex_rxn_ids = String[]
    for rxn_id in reactions(model)
        if any(startswith(rxn_id, x) for x in ex_prefixes) && !any(occursin(x, rxn_id) for x in biomass_strings) # found exchange reaction
            push!(ex_rxn_ids, rxn_id)
            continue
        elseif !exclude_biomass && any([occursin(x, rxn_id) for x in biomass_strings]) # biomass
            push!(ex_rxn_ids, rxn_id)    
        end
    end

    return ex_rxn_ids
end

"""
    find_exchange_metabolites(
        model::StandardModel;
        exclude_biomass = false,
        biomass_strings =  _constants.biomass_strings,
        ex_prefixes = _constants.exchange_prefixes,
    )::Dict{String, Dict{String, Float64}}

Return a dictionary mapping exchange reaction ids to exchange metabolites in a
dictionary. Identifies these reactions based on their prefixes, see
`_constants.exchange_prefixes` for the default list (includes prefixes like
"EX_", etc.). If `exclude_biomass` is true then the biomass reaction is not
returned as in this list, otherwise looks for the biomass reaction by checking
if a reaction contains a string in the list `biomass_strings`, see
`_constants.biomass_strings` for the default list (includes strings like
"biomass", etc.).
"""
function find_exchange_metabolites(
    model::M;
    exclude_biomass = false,
    biomass_strings =  _constants.biomass_strings,
    ex_prefixes = _constants.exchange_prefixes,
)::Dict{String, Dict{String, Float64}} where M <: MetabolicModel
    ex_rxns = find_exchange_reactions(
        model,
        exclude_biomass = exclude_biomass,
        biomass_strings = biomass_strings,
        ex_prefixes = ex_prefixes,
    )
    ex_rxn_met = Dict{String, Dict{String, Float64}}()
    for ex_rxn in ex_rxns
        ex_rxn_met[ex_rxn] = reaction_equation(model, ex_rxn)
    end
    return ex_rxn_met
end
