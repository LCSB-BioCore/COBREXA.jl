"""
    join

Join together models in `models`.
NB: same namespace must be used for both reactions and metabolites.
"""
function Base.join(models::Vector{MetabolicModel}; objective_weights=[])
    n_total_model_reactions = [n_reactions(model) for model in models]
    n_total_model_metabolites = [n_metabolites(model) for model in models]

    metabolic_block = spzeros(sum(n_total_model_metabolites), sum(n_total_model_reactions))
    lbs = spzeros(sum(n_total_model_reactions))
    ubs = spzeros(sum(n_total_model_reactions))

    last_rxn = cumsum(n_total_model_reactions)
    last_met = cumsum(n_total_model_metabolites)
    first_rxn = [1; cumsum(1 .+ n_total_model_reactions)]
    first_met = [1; cumsum(1 .+ n_total_model_metabolites)]
    for (i, fm, lm, fr, lr) in zip(1:length(models),first_met, last_met, first_rxn, last_rxn)
        metabolic_block[fm:lm, fr:lr] .= stoichiometry(model[i])
        tlbs, tubs = bounds(models[i])
        lbs[fr:lr] .= tlbs
        ubs[fr:lr] .= ulbs 
    end

    environment_block = spzeros(n_total_exchanges, sum(n_total_model_reactions))

    if isempty(objective_weights)
        return vcat(metabolic_block, environment_block), lbs, ubs
    else
        objective_column = spzeros(n_total_exchanges+n_total_model_metabolites)
        return hcat(vcat(metabolic_block, environment_block), objective_column), lbs, ubs
    end
end
