"""
    join

Join together models in `models` using the exchange metabolites with ids in `exchange_ids`.
"""
function Base.join(models::Vector{M}, exchange_rxn_ids::Vector{String}, exchange_met_ids::Vector{String}; 
    objective_inds = [], objective_weights=[], species_names=[]
    ) where {M<:MetabolicModel}
    
    n_total_model_reactions = [n_reactions(model) for model in models]
    n_total_model_metabolites = [n_metabolites(model) for model in models]

    metabolic_block = spzeros(sum(n_total_model_metabolites), sum(n_total_model_reactions))
    lbs = spzeros(sum(n_total_model_reactions) + length(exchange_rxn_ids))
    ubs = spzeros(sum(n_total_model_reactions) + length(exchange_rxn_ids))
    metabolic_rxns = Array{String, 1}(undef, sum(n_total_model_reactions))
    metabolic_mets = Array{String, 1}(undef, sum(n_total_model_metabolites))

    environment_block = spzeros(length(exchange_met_ids), sum(n_total_model_reactions))

    last_rxn = cumsum(n_total_model_reactions)
    last_met = cumsum(n_total_model_metabolites)
    first_rxn = [1; cumsum(1 .+ n_total_model_reactions[1:end-1])]
    first_met = [1; cumsum(1 .+ n_total_model_metabolites[1:end-1])]
    for (i, fm, lm, fr, lr) in zip(1:length(models),first_met, last_met, first_rxn, last_rxn)
        species = isempty(species_names) ? "species" : species_names[i]
        metabolic_block[fm:lm, fr:lr] .= stoichiometry(models[i])
        tlbs, tubs = bounds(models[i])
        lbs[fr:lr] .= tlbs
        ubs[fr:lr] .= tubs
        metabolic_mets[fm:lm] .= "$(species)_$(i)_".*metabolites(models[i])
        metabolic_rxns[fr:lr] .= "$(species)_$(i)_".*reactions(models[i])
        exchange_rxn_inds = indexin(exchange_rxn_ids, reactions(models[i]))
        exchange_met_inds = indexin(exchange_met_ids, metabolites(models[i]))

        for (j, ex_met, ex_rxn) in zip(1:length(exchange_rxn_inds), exchange_met_inds, exchange_rxn_inds) # each exchange rxn has one exchange met
            (isnothing(ex_rxn) || isnothing(ex_met)) && continue
            environment_block[j, fr + ex_rxn - 1] = -metabolic_block[fm:lm, fr:lr][ex_met, ex_rxn]
        end
    end

    S = vcat(metabolic_block, environment_block)
    environment_exchange_block = vcat(spzeros(size(S, 1)-length(exchange_met_ids), length(exchange_rxn_ids)), spdiagm(ones(length(exchange_rxn_ids))))
    S = hcat(S, environment_exchange_block)
    all_rxns = [metabolic_rxns; "env_".*exchange_rxn_ids]
    all_mets = [metabolic_mets; "env_".*exchange_met_ids]


    
    return CoreModel(
            S, 
            spzeros(size(S, 1)), 
            spzeros(size(S, 2)), 
            lbs, 
            ubs, 
            all_rxns, 
            all_mets,
            )
end


function all_boundaries(model::M) where {M<:MetabolicModel}
    boundary_mets = String[]
    boundary_rxns = String[]
    S = stoichiometry(model)
    rxns = reactions(model)
    mets = metabolites(model)
    for i in 1:size(S, 2)
        j, b = findnz(S[:, i])
        if length(j) == 1
            push!(boundary_mets, mets[first(j)])
            push!(boundary_rxns, rxns[i])
        end
    end
    return boundary_rxns, boundary_mets
end
