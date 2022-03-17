"""
    looks_like_exchange_reaction(rxn_id::String;
        exclude_biomass = false,
        biomass_strings = _constants.biomass_strings,
        exchange_prefixes = _constants.exchange_prefixes,
    )

A predicate that matches reaction identifiers that look like
exchange or biomass reactions, given the usual naming schemes in common model
repositories. Exchange reactions are identified based on matching prefixes in
the set `exchange_prefixes` and biomass reactions are identified by looking for
occurences of `biomass_strings` in the reaction id.

Also see [`find_exchange_reactions`](@ref).

# Example
```
findall(looks_like_exchange_reaction, reactions(model)) # returns indices
filter(looks_like_exchange_reaction, reactions(model)) # returns Strings

# to use the optional arguments you need to expand the function's arguments
# using an anonymous function
findall(x -> looks_like_exchange_reaction(x; exclude_biomass=true), reactions(model)) # returns indices
filter(x -> looks_like_exchange_reaction(x; exclude_biomass=true), reactions(model)) # returns Strings
```
"""
function looks_like_exchange_reaction(
    rxn_id::String;
    exclude_biomass = false,
    biomass_strings = _constants.biomass_strings,
    exchange_prefixes = _constants.exchange_prefixes,
)::Bool
    any(startswith(rxn_id, x) for x in exchange_prefixes) &&
        !(exclude_biomass && any(occursin(x, rxn_id) for x in biomass_strings))
end

"""
    find_exchange_reactions(m::MetabolicModel; kwargs...)

Shortcut for finding exchange reaction indexes in a model; arguments are
forwarded to [`looks_like_exchange_reaction`](@ref).
"""
find_exchange_reactions(m::MetabolicModel; kwargs...) =
    findall(id -> looks_like_exchange_reaction(id; kwargs...), reactions(m))

"""
    find_exchange_reaction_ids(m::MetabolicModel; kwargs...)

Shortcut for finding exchange reaction identifiers in a model; arguments are
forwarded to [`looks_like_exchange_reaction`](@ref).
"""
find_exchange_reaction_ids(m::MetabolicModel; kwargs...) =
    filter(id -> looks_like_exchange_reaction(id, kwargs...), reactions(m))

"""
    looks_like_biomass_reaction(rxn_id::String;
        exclude_exchanges = false,
        exchange_prefixes = _constants.exchange_prefixes,
        biomass_strings = _constants.biomass_strings,
    )::Bool

A predicate that matches reaction identifiers that look like biomass reactions.
Biomass reactions are identified by looking for occurences of `biomass_strings`
in the reaction id. If `exclude_exchanges` is set, the strings that look like
exchanges (from [`looks_like_exchange_reaction`](@ref)) will not match.

# Example
```
filter(looks_like_biomass_reaction, reactions(model)) # returns strings
findall(looks_like_biomass_reaction, reactions(model)) # returns indices
```
"""
function looks_like_biomass_reaction(
    rxn_id::String;
    exclude_exchanges = false,
    exchange_prefixes = _constants.exchange_prefixes,
    biomass_strings = _constants.biomass_strings,
)::Bool
    any(occursin(x, rxn_id) for x in biomass_strings) &&
        !(exclude_exchanges && any(startswith(rxn_id, x) for x in exchange_prefixes))
end

"""
    find_biomass_reactions(m::MetabolicModel; kwargs...)

Shortcut for finding biomass reaction indexes in a model; arguments are
forwarded to [`looks_like_biomass_reaction`](@ref).
"""
find_biomass_reactions(m::MetabolicModel; kwargs...) =
    findall(id -> looks_like_biomass_reaction(id; kwargs...), reactions(m))

"""
    find_biomass_reaction_ids(m::MetabolicModel; kwargs...)

Shortcut for finding biomass reaction identifiers in a model; arguments are
forwarded to [`looks_like_biomass_reaction`](@ref).
"""
find_biomass_reaction_ids(m::MetabolicModel; kwargs...) =
    filter(id -> looks_like_biomass_reaction(id; kwargs...), reactions(m))

"""
    looks_like_extracellular_metabolite(rxn_id::String;
        extracellular_suffixes = _constants.extracellular_suffixes,
    )::Bool

A predicate that matches metabolite identifiers that look like they are extracellular
metabolites. Extracellular metabolites are identified by `extracellular_suffixes` at the end of the
metabolite id.

# Example
```
filter(looks_like_extracellular_metabolite, metabolites(model)) # returns strings
findall(looks_like_extracellular_metabolite, metabolites(model)) # returns indices
```
"""
function looks_like_extracellular_metabolite(
    met_id::String;
    extracellular_suffixes = _constants.extracellular_suffixes,
)::Bool
    any(endswith(met_id, x) for x in extracellular_suffixes)
end

"""
    find_extracellular_metabolites(m::MetabolicModel; kwargs...)

Shortcut for finding extracellular metabolite indexes in a model; arguments are
forwarded to [`looks_like_extracellular_metabolite`](@ref).
"""
find_extracellular_metabolites(m::MetabolicModel; kwargs...) =
    findall(id -> looks_like_extracellular_metabolite(id; kwargs...), metabolites(m))
"""
    find_extracellular_metabolite_ids(m::MetabolicModel; kwargs...)

Shortcut for finding extracellular metabolite identifiers in a model; arguments are
forwarded to [`looks_like_extracellular_metabolite`](@ref).
"""
find_extracellular_metabolite_ids(m::MetabolicModel; kwargs...) =
    findall(id -> looks_like_extracellular_metabolite(id; kwargs...), metabolites(m))
