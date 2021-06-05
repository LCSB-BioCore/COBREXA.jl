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
Returns indices of exchange reactions.
Exchange reactions are identified based on most commonly used prefixes.
"""
function find_exchange_reactions(
    model::CoreModel;
    exclude_biomass = false,
    biomass_str::String = "biomass",
    exc_prefs = ["EX_"; "Exch_"; "Ex_"; "R_EX_"],
)
    is_exc = falses(n_reactions(model))
    for pref in exc_prefs
        is_exc = is_exc .| startswith.(model.rxns, pref)
    end
    exc_inds = findall(is_exc)
    if exclude_biomass
        biom_inds = findall(x -> occursin(biomass_str, x), model.rxns)
        exc_inds = setdiff(exc_inds, biom_inds)
    end
    return exc_inds
end

"""
Returns indices of exchanged metabolites, ie, the outermost metabolites in the network
In practice returns the metabolites consumed by the reactions given by `find_exchange_reactions`
and if called with the same arguments, the two outputs correspond.
"""
function find_exchange_metabolites(
    model::CoreModel;
    exclude_biomass = false,
    biomass_str::String = "biomass",
    exc_prefs = ["EX_"; "Exch_"; "Ex_"; "R_EX_"],
)
    exc_rxn_inds = find_exchange_reactions(
        model,
        exclude_biomass = exclude_biomass,
        biomass_str = biomass_str,
        exc_prefs = exc_prefs,
    )
    exc_met_inds = [findfirst(x -> x == -1, model.S[:, j]) for j in exc_rxn_inds]
    return exc_met_inds
end
