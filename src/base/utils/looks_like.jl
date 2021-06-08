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
