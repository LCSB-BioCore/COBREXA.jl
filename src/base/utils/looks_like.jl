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
    filter(id -> looks_like_biomass_reaction(id, kwargs...), reactions(m))

"""
    looks_like_internal_reaction(
        rxn_id::String;
        exchange_prefixes = _constants.exchange_prefixes,
        biomass_strings = _constants.biomass_strings,
    )::Bool

A predicate that matches reaction identifiers that look like
internal reactions, i.e. reactions that are neither exchange nor biomass reactions.
Exchange reactions are identified based on matching prefixes in
the set `exchange_prefixes` and biomass reactions are identified by looking for
occurences of `biomass_strings` in the reaction id.

Also see [`find_internal_reactions`](@ref).

# Example
```
findall(looks_like_internal_reaction, reactions(model)) # returns indices
filter(looks_like_internal_reaction, reactions(model)) # returns Strings

# to use the optional arguments you need to expand the function's arguments
# using an anonymous function
findall(x -> looks_like_internal_reaction(x; exchange_prefixes=["EX_", "R_EX_"]), reactions(model)) # returns indices
filter(x -> looks_like_internal_reaction(x; exchange_prefixes=["EX_", "R_EX_"]), reactions(model)) # returns Strings
```
"""
function looks_like_internal_reaction(
    rxn_id::String;
    exchange_prefixes = _constants.exchange_prefixes,
    biomass_strings = _constants.biomass_strings,
)::Bool
    !looks_like_biomass_reaction(
        rxn_id;
        exchange_prefixes = exchange_prefixes,
        biomass_strings = biomass_strings,
    ) &&
        !looks_like_exchange_reaction(
            rxn_id;
            exchange_prefixes = exchange_prefixes,
            biomass_strings = biomass_strings,
        )
end

"""
    find_internal_reactions(m::MetabolicModel; kwargs...)

Shortcut for finding internal reaction indices in a model; arguments are
forwarded to [`looks_like_internal_reaction`](@ref).
"""
find_internal_reactions(m::MetabolicModel; kwargs...) =
    findall(id -> looks_like_internal_reaction(id; kwargs...), reactions(m))

"""
    find_internal_reaction_ids(m::MetabolicModel; kwargs...)

Shortcut for finding internal reaction identifiers in a model; arguments are
forwarded to [`looks_like_internal_reaction`](@ref).
"""
find_internal_reaction_ids(m::MetabolicModel; kwargs...) =
    filter(id -> looks_like_internal_reaction(id, kwargs...), reactions(m))

"""
    looks_like_exchange_metabolite(rxn_id::String;
        exchange_suffixes = _constants.exchange_suffixes,
        )::Bool

A predicate that matches metabolite identifiers that look like involved in
exchange reactions. Exchange metabolites are identified by `exchange_suffixes`
at the end of the metabolite id.

# Example
```
filter(looks_like_exchange_metabolite, metabolites(model)) # returns strings
findall(looks_like_exchange_metabolite, metabolites(model)) # returns indices
```
"""
function looks_like_exchange_metabolite(
    met_id::String;
    exchange_suffixes = _constants.exchange_suffixes,
)::Bool
    any(endswith(met_id, x) for x in exchange_suffixes)
end

"""
    find_exchange_metabolites(m::MetabolicModel; kwargs...)

Shortcut for finding exchange metabolite indexes in a model; arguments are
forwarded to [`looks_like_exchange_metabolite`](@ref).
"""
find_exchange_metabolites(m::MetabolicModel; kwargs...) =
    findall(id -> looks_like_exchange_metabolite(id; kwargs...), metabolites(m))
"""
    find_exchange_metabolite_ids(m::MetabolicModel; kwargs...)

Shortcut for finding exchange metabolite identifiers in a model; arguments are
forwarded to [`looks_like_exchange_metabolite`](@ref).
"""
find_exchange_metabolite_ids(m::MetabolicModel; kwargs...) =
    findall(id -> looks_like_exchange_metabolite(id; kwargs...), metabolites(m))
