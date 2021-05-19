"""
add_objective!(cmodel, objective_mets, objective_weights, objective_column)


"""
function add_objective!()

end

"""
append!()

Append `model` to `cmodel` where `cmodel` is a pre-existing community model with `exchange_rxn_ids` and
`exchange_met_ids`. If an objective function has already been assigned then supply its column index in `objective_col`
and the metabolites used by the objective in `objective_rows` as well as the weight to assign the new
"""
function Base.append!(
    cmodel::CoreModel,
    model::M,
    exchange_rxn_ids::Vector{String},
    exchange_met_ids::Vector{String};
    species_name = "",
    objective_ind = 0,
) where {M<:MetabolicModel}
    cmI, cmJ, cmV = findnz(cmodel)
    mI, mJ, mV = findnz(model)

    nI = zeros(Int, length(cmI) + length(mI))
    nJ = zeros(Int, length(cmJ) + length(mJ))
    nV = zeros(Int, length(cmV) + length(mV))

    # find location of environmental exchange reactions
    ex_mets = indexin(exchange_met_ids, metabolites(cmodel))
    ex_rxns = indexin(exchange_rxn_ids, reactions(cmodel))


end

"""
    Base.join(models::Vector{M}, 
        exchange_rxn_ids::Vector{String}, 
        exchange_met_ids::Vector{String}; 
        add_biomass_objective=false, 
        biomass_ids::Vector{String}, 
        species_names=String[]
    ) 

Return a `CoreModel` representing the community model of `models` joined
through their `exchange_rxn_ids` and `exchange_met_ids`. These exchange
reactions and metabolites are converted into environmental metabolites that link the models.
Optionally specify `species_names` to append a specific name to each reaction
and metabolite of an organism for easier reference. Note, the bounds of the 
environmental variables are all set to zero. Thus, to run a simulation
you need to constrain them appropriately. All the other bounds are inherited
from the models used to construct the community model.

If `add_biomass_objective` is true then `biomass_ids` needs to be supplied as
well. This creates a model with an extra reaction added to the end of the
stoichiometric matrix (last column) that can be assigned as the objective
reaction. It also creates biomass "metabolites" that can be used in this objective
reaction. Note, this reaction is unspecified, further action needs to be taken
to specify it, e.g. assign weights to the last column of the stoichiometric matrix
in the rows corresponding to the biomass metabolites.
"""
function Base.join(
    models::Vector{M},
    exchange_rxn_ids::Vector{String},
    exchange_met_ids::Vector{String};
    add_biomass_objective = true,
    biomass_ids = String[],
    species_names = String[],
) where {M<:MetabolicModel}

    if add_biomass_objective && isempty(biomass_ids)
        throw(
            DomainError(
                "biomass_ids",
                "Please add supply the string ids of the biomass functions when `add_biomass_objective` is true.",
            ),
        )
    end

    # get offsets to construct community S
    reaction_lengths = [n_reactions(model) for model in models]
    metabolite_lengths = [n_metabolites(model) for model in models]
    reaction_offset = [0; cumsum(reaction_lengths[1:end-1])]
    metabolite_offset = [0; cumsum(metabolite_lengths[1:end-1])]

    # get each model's S matrix (needed for the size calculations)
    stoichs = [stoichiometry(model) for model in models]
    nnzs = [findnz(stoich) for stoich in stoichs] # nonzero indices per model

    # size calculations
    column_add = add_biomass_objective ? 1 : 0 # objective rxn
    row_add = add_biomass_objective ? length(models) : 0 # biomass as metabolites
    nnz_add = add_biomass_objective ? (1 + length(models)) : 0
    nnz_total =
        sum(length(first(nnz)) for nnz in nnzs) +
        length(models) * length(exchange_rxn_ids) +
        length(exchange_met_ids) + nnz_add
    n_reactions_metabolic = sum(reaction_lengths)
    n_reactions_total = n_reactions_metabolic + length(exchange_rxn_ids) + column_add
    n_metabolites_metabolic = sum(metabolite_lengths)
    n_metabolites_total = n_metabolites_metabolic + length(exchange_met_ids) + row_add

    # Create empty storage vectors
    lbs = spzeros(n_reactions_total)
    ubs = spzeros(n_reactions_total)
    rxns = Array{String,1}(undef, n_reactions_total)
    mets = Array{String,1}(undef, n_metabolites_total)
    I = Array{Int,1}(undef, nnz_total)
    J = Array{Int,1}(undef, nnz_total)
    V = Array{Float64,1}(undef, nnz_total)

    # build metabolic components
    kstart = 1
    for i = 1:length(models)
        kend = kstart + length(nnzs[i][3]) - 1
        rng = kstart:kend
        I[rng] .= nnzs[i][1] .+ metabolite_offset[i]
        J[rng] .= nnzs[i][2] .+ reaction_offset[i]
        V[rng] .= nnzs[i][3]
        kstart += length(nnzs[i][3])
    end
    # build environmental components
    for i = 1:length(models)
        exchange_rxn_inds = indexin(exchange_rxn_ids, reactions(models[i]))
        for (n, ex_rxn) in enumerate(exchange_rxn_inds) # each exchange rxn has one exchange met
            isnothing(ex_rxn) && continue
            # connect environmental metabolite with exchange metabolite
            I[kstart] = n_metabolites_metabolic + n
            J[kstart] = ex_rxn + reaction_offset[i]
            V[kstart] = nnzs[i][3][ex_rxn]
            kstart += 1
        end
    end
    # build diagonal environmental exchanges
    for i = 1:length(exchange_rxn_ids)
        I[kstart] = n_metabolites_metabolic + i
        J[kstart] = n_reactions_metabolic + i
        V[kstart] = -1.0
        kstart += 1
    end

    if add_biomass_objective
        n_before_biomass_row = n_metabolites_metabolic + length(exchange_met_ids)
        for i = 1:length(models)
            biomass_ind = first(indexin([biomass_ids[i]], reactions(models[i])))
            I[kstart] = i + n_before_biomass_row
            J[kstart] = biomass_ind + reaction_offset[i]
            V[kstart] = 1.0
            kstart += 1
        end
    end
    S = sparse(I[1:kstart-1], J[1:kstart-1], V[1:kstart-1], n_metabolites_total, n_reactions_total) # could be that some microbes don't have all the exchanges

    reaction_cumsum = cumsum(reaction_lengths)
    metabolite_cumsum = cumsum(metabolite_lengths)
    for i = 1:length(models)
        species = isempty(species_names) ? "species_$(i)" : species_names[i]
        tlbs, tubs = bounds(models[i])
        lbs[reaction_offset[i]+1:reaction_cumsum[i]] .= tlbs
        ubs[reaction_offset[i]+1:reaction_cumsum[i]] .= tubs
        rxns[reaction_offset[i]+1:reaction_cumsum[i]] =
            "$(species)_" .* reactions(models[i])
        mets[metabolite_offset[i]+1:metabolite_cumsum[i]] =
            "$(species)_" .* metabolites(models[i])
    end
    mets[metabolite_cumsum[end]+1:metabolite_cumsum[end]+length(exchange_met_ids)] .= exchange_met_ids
    rxns[reaction_cumsum[end]+1:reaction_cumsum[end]+length(exchange_rxn_ids)] .= exchange_rxn_ids
    
    if add_biomass_objective
        rxns[end] = "community_biomass"
        for i=1:length(models)
            species = isempty(species_names) ? "species_$(i)" : species_names[i]
            mets[end-length(biomass_ids)+i] = "$(species)_".*biomass_ids[i]
        end
    end

    return CoreModel(S, spzeros(size(S, 1)), spzeros(size(S, 2)), lbs, ubs, rxns, mets)
end

"""
    all_boundaries(model::M) where {M<:MetabolicModel}

Return `boundary_rxns` and `boundary_mets` given `model`. Boundary reactions
are reactions with a single stoichiometric coefficient in a column in the
stoichiometric matrix, and boundary metabolites are the corresponding row/metabolite
for that column.
"""
function all_boundaries(model::M) where {M<:MetabolicModel}
    boundary_mets = String[]
    boundary_rxns = String[]
    S = stoichiometry(model)
    rxns = reactions(model)
    mets = metabolites(model)
    for i = 1:size(S, 2)
        j, b = findnz(S[:, i])
        if length(j) == 1
            push!(boundary_mets, mets[first(j)])
            push!(boundary_rxns, rxns[i])
        end
    end
    return boundary_rxns, boundary_mets
end
