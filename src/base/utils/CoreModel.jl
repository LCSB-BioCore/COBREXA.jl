Base.isequal(model1::CoreModel, model2::CoreModel) =
    isequal(model1.S, model2.S) &&
    isequal(model1.b, model2.b) &&
    isequal(model1.c, model2.c) &&
    isequal(model1.xl, model2.xl) &&
    isequal(model1.xu, model2.xu) &&
    isequal(model1.rxns, model2.rxns) &&
    isequal(model1.mets, model2.mets)

Base.copy(model::CoreModel) =
    CoreModel(model.S, model.b, model.c, model.xl, model.xu, model.rxns, model.mets)

Base.isequal(model1::CoreModelCoupled, model2::CoreModelCoupled) =
    isequal(model1.lm, model2.lm) &&
    isequal(model1.C, model2.C) &&
    isequal(model1.cl, model2.cl) &&
    isequal(model1.cu, model2.cu)

Base.copy(model::CoreModelCoupled) = CoreModelCoupled(model.lm, model.C, model.cl, model.cu)

"""
    find_exchange_reactions(
        model::CoreModel;
        exclude_biomass = false,
        biomass_strings = _constants.biomass_strings,
        ex_prefixes = _constants.exchange_prefixes,
    )::Vector{String}

Returns indices of exchange reactions. Exchange reactions are identified based
on most commonly used prefixes, these prefixes can be set using `ex_prefixes`.
See `_constants.biomass_strings` for the default list. The biomass reaction is
also added if it can be identified by looking for the strings listed in
`biomass_strings` in the reaction ids. If `exclude_biomass` is true then this
does not occur.

Note: biomass exchange reactions are counted as exchange reactions and will NOT
be excluded when `exclude_biomass` is true.
"""
function find_exchange_reactions(
    model::CoreModel;
    exclude_biomass = false,
    biomass_strings = _constants.biomass_strings,
    ex_prefixes = _constants.exchange_prefixes,
)::Vector{Int}
    ex_inds = Int[]
    for (i, rxn_id) in enumerate(reactions(model))
        if any([startswith(rxn_id, x) for x in ex_prefixes]) # exchange reactions
            push!(ex_inds, i)
            continue
        elseif !exclude_biomass && any([occursin(x, rxn_id) for x in biomass_strings]) # biomass
            push!(ex_inds, i)
        end
    end
    return ex_inds
end

"""
    find_exchange_metabolites(
        model::CoreModel;
        exclude_biomass = false,
        biomass_strings = _constants.biomass_strings,
        ex_prefixes = _constants.exchange_prefixes,
    )::Vector{String}

Returns a dictionary mapping indices of exchange reactions to dictionaries of
exchange metabolites where the metabolite dictionary corresponds to metabolite
indices mapped to stoichiometric coefficients. Exchange reactions are identified
based on most commonly used prefixes, these prefixes can be set using
`ex_prefixes`. See `_constants.biomass_strings` for the default list. The
biomass reaction is also added if it can be identified by looking for the
strings listed in `biomass_strings` in the reaction ids. If `exclude_biomass` is
true then this does not occur.

Note: biomass exchange reactions are counted as exchange reactions and will NOT
be excluded when `exclude_biomass` is true.
"""
function find_exchange_metabolites(
    model::CoreModel;
    exclude_biomass = false,
    biomass_strings = _constants.biomass_strings,
    ex_prefixes = _constants.exchange_prefixes,
)::Dict{Int, Dict{Int, Float64}}
    exc_rxn_inds = find_exchange_reactions(
        model,
        exclude_biomass = exclude_biomass,
        biomass_strings = biomass_strings,
        ex_prefixes = ex_prefixes,
    )
    exc_met_inds = Dict(exc_rxn_ind =>
        Dict(k=>v for (k, v) in zip(findnz(model.S[:, exc_rxn_ind])...)) for exc_rxn_ind in exc_rxn_inds
    )
    return exc_met_inds
end
