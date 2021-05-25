"""
add_objective!(cmodel, objective_mets, objective_weights, objective_column)


"""
function add_objective!()

end



"""
push!()

Append `model` to `cmodel` where `cmodel` is a pre-existing community model with `exchange_rxn_ids` and
`exchange_met_ids`. If an objective function has already been assigned then supply its column index in `objective_col`
and the metabolites used by the objective in `objective_rows` as well as the weight to assign the new
"""
function Base.push!(
    cmodel::CoreModel,
    model::M,
    exchange_met_ids::Vector{String},
    exchange_rxn_ids::Vector{String},
    species_name,
    has_biomass_objective;
    biomass_id = "",
    ) where {M<:MetabolicModel}

    if has_biomass_objective && biomass_id == ""
        throw(
            DomainError(
                "Argument required.",
                "The community uses a biomass objective function, please supply the objective id of the model you want to add.",
            ),
        )
    end

    n_cmodel_rows, n_cmodel_cols = size(stoichiometry(cmodel))
    Iadd, Jadd, Vadd = findnz(stoichiometry(model))
    add_row = has_biomass_objective ? 1 : 0
    n_metabolites_total =
        n_reactions_total =
        # shift to fit into bigger model
            Iadd .+= n_cmodel_rows
    Jadd .+= n_cmodel_cols

    exchange_met_community_inds = indexin(exchange_met_ids, metabolites(cmodel))
    if any(isnothing.(exchange_met_community_inds))
        throw(
            DomainError(
                "Exchange metabolite not found.",
                "Exchange metabolite not found in community model.",
            ),
        )
    end
    exchange_rxn_model_inds = indexin(exchange_rxn_ids, reactions(model))

    # when adding a single model not that many reactions, push! okay?
    for i = 1:length(exchange_met_community_inds)
        isnothing(exchange_rxn_model_inds[i]) && continue
        push!(Iadd, n_cmodel_rows + exchange_met_community_inds[i])
        push!(Jadd, n_cmodel_cols + exchange_rxn_model_inds[i])
        push!(Vadd, 1.0)
    end

    # if has_biomass_objective
    #     biomass_ind = first(indexin([biomass_id]), reactions(model))

    #     push!(Iadd, n_cmodel_rows + )
    #     push!(Jadd, n_cmodel_cols + )
    #     push!(Vadd, 1.0)
    # end

    I, J, V = findnz(stoichiometry(cmodel))
    I = [I; Iadd]
    J = [J; Jadd]
    V = [V; Vadd]
    S = sparse(I, J, V, n_metabolites_total, n_reactions_total)

    lbsadd, ubsadd = bounds(model)
    lbs = [lbs; lbsadd]
    ubs = [lbs; ubsadd]

    rxnsadd = "$(species_name)_".reactions(model)
    if has_biomass_objective
        metsadd =
            ["$(species_name)_" .* metabolites(model); "$(species_name)_" * biomass_id]
    else
        metsadd = "$(species_name)_".metabolites(model)
    end
    rxns = [rxns; rxnsadd]
    mets = [mets; metsadd]

    return CoreModel(S, spzeros(size(S, 1)), spzeros(size(S, 2)), lbs, ubs, rxns, mets)
end

"""
    join(models::Vector{M}, 
        exchange_rxn_ids::Vector{String}, 
        exchange_met_ids::Vector{String}; 
        add_biomass_objective=false, 
        biomass_ids::Vector{String}, 
        species_names=String[]
    ) 

Return a `CoreModel` representing the community model of `models` joined through
their `exchange_rxn_ids` and `exchange_met_ids`. These exchange reactions and
metabolites link to environmental metabolites and reactions.
Optionally specify `species_names` to append a specific name to each reaction
and metabolite of an organism for easier reference (default is `species_i` for
each model index i in `models`). Note, the bounds of the environmental variables
are all set to zero. Thus, to run a simulation you need to constrain them
appropriately. All the other bounds are inherited from the models used to
construct the community model.

If `add_biomass_objective` is true then `biomass_ids` needs to be supplied as
well. This creates a model with an extra reaction added to the end of the
stoichiometric matrix (last column) that can be assigned as the objective
reaction. It also creates biomass "metabolites" that can be used in this
objective reaction. Note, this reaction is unspecified, further action needs to
be taken to specify it, e.g. assign weights to the last column of the
stoichiometric matrix in the rows corresponding to the biomass metabolites.

To further clarify how this `join` works. Suppose you have 2 organisms with
stoichiometric matrices S₁ and S₂ and you want to link them with
`exchange_met_ids = [em₁, em₂, em₃, ...]` and `exchange_rxn_ids = [er₁, er₂,
er₃, ...]`. Then a new community stoichiometric matrix is constructed that looks
like this:
```
            _      er₁  er₂  er₃  ...  b_
           |S₁                           |
           |   S₂                        |
        em₁|                             |
S   =   em₂|                             |
        em₃|                             |
        ...|                             |
        bm₁|                             |  
        bm₂|_                           _|

```
The exchange reactions in each model get linked to environmental metabolites,
`emᵢ`, and these get linked to environmental exchanges, `erᵢ`. These `erᵢ`
behave like normal single organism exchange reactions. When
`add_biomass_objective` is true each model's biomass becomes a pseudo-metabolite
(`bmᵢ`). These can be weighted in column `b`, called the `community_biomass`
reaction in the community model, if desired. Refer to the tutorial if this is
unclear.

# Example
```
model_path = joinpath("..","models","e_coli_core.json")
m1 = load_model(StandardModel, model_path)
model_path = joinpath("iML1515.xml")
m2 = load_model(StandardModel, model_path)

exchange_rxn_ids, exchange_met_ids = all_boundaries(m2)
biomass_ids = ["BIOMASS_Ecoli_core_w_GAM","R_BIOMASS_Ec_iML1515_core_75p37M"]
community = COBREXA.join([m1, m2], exchange_rxn_ids, exchange_met_ids; add_biomass_objective=true, biomass_ids=biomass_ids, species_names=["Core", "iML1515"])
```
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

    # build metabolic components, block diagonals
    kstart = 1
    for i = 1:length(models)
        kend = kstart + length(nnzs[i][3]) - 1
        rng = kstart:kend
        I[rng] .= nnzs[i][1] .+ metabolite_offset[i]
        J[rng] .= nnzs[i][2] .+ reaction_offset[i]
        V[rng] .= nnzs[i][3]
        kstart += length(nnzs[i][3])
    end

    # build environmental - exchange links
    for i = 1:length(models)
        exchange_rxn_inds = indexin(exchange_rxn_ids, reactions(models[i]))
        exchange_met_inds = indexin(exchange_met_ids, metabolites(models[i]))
        for (n, (ex_rxn, ex_met)) in enumerate(zip(exchange_rxn_inds, exchange_met_inds)) # each exchange rxn has one exchange met
            isnothing(ex_rxn) && continue
            isnothing(ex_met) && continue
            # connect environmental metabolite with exchange metabolite
            I[kstart] = n_metabolites_metabolic + n
            J[kstart] = ex_rxn + reaction_offset[i]
            V[kstart] = -stoichs[i][ex_met, ex_rxn] # ex is normally negative, make positive
            kstart += 1
        end
    end

    # # build diagonal environmental exchanges
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

    S = sparse(
        I[1:kstart-1],
        J[1:kstart-1],
        V[1:kstart-1],
        n_metabolites_total,
        n_reactions_total,
    ) # could be that some microbes don't have all the exchanges, hence kstart-1

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
    mets[metabolite_cumsum[end]+1:metabolite_cumsum[end]+length(exchange_met_ids)] .=
        exchange_met_ids
    rxns[reaction_cumsum[end]+1:reaction_cumsum[end]+length(exchange_rxn_ids)] .=
        exchange_rxn_ids

    if add_biomass_objective
        rxns[end] = "community_biomass"
        for i = 1:length(models)
            species = isempty(species_names) ? "species_$(i)" : species_names[i]
            mets[end-length(biomass_ids)+i] = "$(species)_" .* biomass_ids[i]
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
