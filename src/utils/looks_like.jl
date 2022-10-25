"""
$(TYPEDSIGNATURES)

A predicate that matches reaction identifiers that look like
exchange or biomass reactions, given the usual naming schemes in common model
repositories. Exchange reactions are identified based on matching prefixes in
the set `exchange_prefixes` and biomass reactions are identified by looking for
occurences of `biomass_strings` in the reaction id.

Also see [`find_exchange_reactions`](@ref).

# Note
While `looks_like_exchange_reaction` is useful for heuristically finding a
reaction, it is preferable to use standardized terms for finding reactions (e.g.
SBO terms). See [`is_exchange_reaction`](@ref) for a more systematic
alternative.

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
    biomass_strings = constants.biomass_strings,
    exchange_prefixes = constants.exchange_prefixes,
)::Bool
    any(startswith(rxn_id, x) for x in exchange_prefixes) &&
        !(exclude_biomass && any(occursin(x, rxn_id) for x in biomass_strings))
end

"""
$(TYPEDSIGNATURES)

Shortcut for finding exchange reaction indexes in a model; arguments are
forwarded to [`looks_like_exchange_reaction`](@ref).
"""
find_exchange_reactions(m::AbstractMetabolicModel; kwargs...) =
    findall(id -> looks_like_exchange_reaction(id; kwargs...), reactions(m))

"""
$(TYPEDSIGNATURES)

Shortcut for finding exchange reaction identifiers in a model; arguments are
forwarded to [`looks_like_exchange_reaction`](@ref).
"""
find_exchange_reaction_ids(m::AbstractMetabolicModel; kwargs...) =
    filter(id -> looks_like_exchange_reaction(id, kwargs...), reactions(m))

"""
$(TYPEDSIGNATURES)

A predicate that matches reaction identifiers that look like biomass reactions.
Biomass reactions are identified by looking for occurences of `biomass_strings`
in the reaction id. If `exclude_exchanges` is set, the strings that look like
exchanges (from [`looks_like_exchange_reaction`](@ref)) will not match.

# Note
While `looks_like_biomass_reaction` is useful for heuristically finding a
reaction, it is preferable to use standardized terms for finding reactions (e.g.
SBO terms). See [`is_biomass_reaction`](@ref) for a more systematic
alternative.

# Example
```
filter(looks_like_biomass_reaction, reactions(model)) # returns strings
findall(looks_like_biomass_reaction, reactions(model)) # returns indices
```
"""
function looks_like_biomass_reaction(
    rxn_id::String;
    exclude_exchanges = false,
    exchange_prefixes = constants.exchange_prefixes,
    biomass_strings = constants.biomass_strings,
)::Bool
    any(occursin(x, rxn_id) for x in biomass_strings) &&
        !(exclude_exchanges && any(startswith(rxn_id, x) for x in exchange_prefixes))
end

"""
$(TYPEDSIGNATURES)

Shortcut for finding biomass reaction indexes in a model; arguments are
forwarded to [`looks_like_biomass_reaction`](@ref).
"""
find_biomass_reactions(m::AbstractMetabolicModel; kwargs...) =
    findall(id -> looks_like_biomass_reaction(id; kwargs...), reactions(m))

"""
$(TYPEDSIGNATURES)

Shortcut for finding biomass reaction identifiers in a model; arguments are
forwarded to [`looks_like_biomass_reaction`](@ref).
"""
find_biomass_reaction_ids(m::AbstractMetabolicModel; kwargs...) =
    filter(id -> looks_like_biomass_reaction(id; kwargs...), reactions(m))

"""
$(TYPEDSIGNATURES)

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
    extracellular_suffixes = constants.extracellular_suffixes,
)::Bool
    any(endswith(met_id, x) for x in extracellular_suffixes)
end

"""
$(TYPEDSIGNATURES)

Shortcut for finding extracellular metabolite indexes in a model; arguments are
forwarded to [`looks_like_extracellular_metabolite`](@ref).
"""
find_extracellular_metabolites(m::AbstractMetabolicModel; kwargs...) =
    findall(id -> looks_like_extracellular_metabolite(id; kwargs...), metabolites(m))

"""
$(TYPEDSIGNATURES)

Shortcut for finding extracellular metabolite identifiers in a model; arguments are
forwarded to [`looks_like_extracellular_metabolite`](@ref).
"""
find_extracellular_metabolite_ids(m::AbstractMetabolicModel; kwargs...) =
    findall(id -> looks_like_extracellular_metabolite(id; kwargs...), metabolites(m))

@_is_reaction_fn "exchange" Identifiers.EXCHANGE_REACTIONS
@_is_reaction_fn "transport" Identifiers.TRANSPORT_REACTIONS
@_is_reaction_fn "biomass" Identifiers.BIOMASS_REACTIONS
@_is_reaction_fn "atp_maintenance" Identifiers.ATP_MAINTENANCE_REACTIONS
@_is_reaction_fn "pseudo" Identifiers.PSEUDOREACTIONS
@_is_reaction_fn "metabolic" Identifiers.METABOLIC_REACTIONS
@_is_reaction_fn "spontaneous" Identifiers.SPONTANEOUS_REACTIONS
