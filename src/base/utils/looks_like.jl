"""
    looks_like_exchange_reaction(rxn_id::String;
        exclude_biomass = false,
        biomass_strings = _constants.biomass_strings,
        exchange_prefixes = _constants.exchange_prefixes,
    )

A predicate function that can be used to identify reactions that are probably 
exchange or biomass reactions. Exchange reactions are identified based on matching prefixes 
in the set `exchange_prefixes` and biomass reactions are identified by looking for occurences 
of `biomass_strings` in the reaction id.

Note: the order of the probable exchange reactions identified here does not necessarily match 
that of the exchange metabolites identified in [`looks_like_exchange_metabolite`](@ref).

# Example
```
filter(looks_like_exchange_reaction, reactions(model)) # returns strings
findall(looks_like_exchange_reaction, reactions(model)) # returns indices

# to use the optional arguments you need to expand the function's arguments using an anonymous function
filter(x -> looks_like_exchange_reaction(x; exclude_biomass=true), reactions(model)) # returns strings
findall(x -> looks_like_exchange_reaction(x; exclude_biomass=true), reactions(model)) # returns indices
```
"""
function looks_like_exchange_reaction(rxn_id::String;
    exclude_biomass = false,
    biomass_strings = _constants.biomass_strings,
    exchange_prefixes = _constants.exchange_prefixes,
)::Bool
    any(startswith(rxn_id, x) for x in exchange_prefixes) && !(exclude_biomass && any(occursin(x, rxn_id) for x in biomass_strings))
end

"""
    looks_like_biomass_reaction(rxn_id::String;
        exclude_exchanges = false,
        exchange_prefixes = _constants.exchange_prefixes,
        biomass_strings = _constants.biomass_strings,
    )::Bool

A predicate function that can be used to identify reactions that are probably 
biomass reactions. Biomass reactions are identified by looking for occurences 
of `biomass_strings` in the reaction id. Can also `exclude_exchanges` that will 
remove occurrences of biomass reactions that have a prefix in `exchange_prefixes`.

# Example
```
filter(looks_like_biomass_reaction, reactions(model)) # returns strings
findall(looks_like_biomass_reaction, reactions(model)) # returns indices
```
"""
function looks_like_biomass_reaction(rxn_id::String;
    exclude_exchanges = false,
    exchange_prefixes = _constants.exchange_prefixes,
    biomass_strings = _constants.biomass_strings,
)::Bool
    any(occursin(x, rxn_id) for x in biomass_strings) && !(exclude_exchanges && any(startswith(rxn_id, x) for x in exchange_prefixes))
end

"""
    looks_like_exchange_metabolite(rxn_id::String;
        exchange_suffixes = _constants.exchange_suffixes,
        )::Bool

A predicate function that can be used to identify metabolites that are probably 
involved in exchange reactions. Exchange metabolites are identified by looking for occurences 
of `exchange_suffixes` at the end of the metabolite id.

Note: the order of the probable exchange metabolites identified here does not necessarily match 
that of the exchange reactions identified in [`looks_like_exchange_reaction`](@ref).

# Example
```
filter(looks_like_exchange_metabolite, metabolites(model)) # returns strings
findall(looks_like_exchange_metabolite, metabolites(model)) # returns indices
```
"""
function looks_like_exchange_metabolite(met_id::String;
    exchange_suffixes = _constants.exchange_suffixes,
)::Bool
    any(endswith(met_id, x) for x in exchange_suffixes)
end
