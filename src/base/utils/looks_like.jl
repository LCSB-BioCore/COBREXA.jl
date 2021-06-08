"""
    looks_like_exchange_reaction(rxn_id::String;
        exclude_biomass = false,
        biomass_strings = _constants.biomass_strings,
        ex_prefixes = _constants.exchange_prefixes,
    )

A predicate function that can be used to identify reactions that are probably 
exchange or biomass reactions. Exchange reactions are identified based on matching prefixes 
in the set `ex_prefixes` and biomass reactions are identified by looking for occurences 
of `biomass_strings` in the reaction id.

Note: the order of the probable exchange reactions identified here does not necessarily match 
that of the exchange metabolites identified in [`looks_like_exchange_metabolite`](@ref).

# Example
```
filter(looks_like_exchange_reaction, reactions(model)) # returns strings
findall(looks_like_exchange_reaction, reactions(model)) # returns indices

# to use the optional arguments you need to instantiate the function
filter(x -> looks_like_exchange_reaction(x; exclude_biomass=true), reactions(model)) # returns strings
findall(x -> looks_like_exchange_reaction(x; exclude_biomass=true), reactions(model)) # returns indices
```
"""
function looks_like_exchange_reaction(rxn_id::String;
    exclude_biomass = false,
    biomass_strings = _constants.biomass_strings,
    ex_prefixes = _constants.exchange_prefixes,
)::Bool
    if exclude_biomass
        if any(startswith(rxn_id, x) for x in ex_prefixes) && !any([occursin(x, rxn_id) for x in biomass_strings])
            return true
        end
    else # don't exclude biomass
        if any(startswith(rxn_id, x) for x in ex_prefixes)
            return true
        end
    end
    return false
end

"""
    looks_like_biomass_reaction(rxn_id::String;
        biomass_strings = _constants.biomass_strings,
    )::Bool

A predicate function that can be used to identify reactions that are probably 
biomass reactions. Biomass reactions are identified by looking for occurences 
of `biomass_strings` in the reaction id.

# Example
```
filter(looks_like_biomass_reaction, reactions(model)) # returns strings
findall(looks_like_biomass_reaction, reactions(model)) # returns indices
```
"""
function looks_like_biomass_reaction(rxn_id::String;
    biomass_strings = _constants.biomass_strings,
)::Bool
    any(occursin(x, rxn_id) for x in biomass_strings) && return true
    return false
end

"""
    looks_like_exchange_metabolite(rxn_id::String;
        ex_suffixes = _constants.exchange_suffixes,
        )::Bool

A predicate function that can be used to identify metabolites that are probably 
involved in exchange reactions. Exchange metabolites are identified by looking for occurences 
of `ex_suffixes` at the end of the metabolite id.

Note: the order of the probable exchange metabolites identified here does not necessarily match 
that of the exchange reactions identified in [`looks_like_exchange_reaction`](@ref).

# Example
```
filter(looks_like_exchange_metabolite, metabolites(model)) # returns strings
findall(looks_like_exchange_metabolite, metabolites(model)) # returns indices
```
"""
function looks_like_exchange_metabolite(met_id::String;
    ex_suffixes = _constants.exchange_suffixes,
)::Bool
    any(endswith(met_id, x) for x in ex_suffixes) && return true
    return false
end
