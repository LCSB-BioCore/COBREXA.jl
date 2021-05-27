"""
add_objective!(community::CoreModel, objective_mets::Vector{String}; objective_weights=Float64[], objective_column_index=0)

Add an objective to the `community` model. Supply the string names of the
objective metabolites in `objective_mets`. Optionally specify the weight to
assign each metabolite in the objective function, if unassigned then equal
weight is assumed. Also, optionally specify whether the objective already exists
in the model by assigning `objective_column_index`. If unassigned then an
objective column will be added, otherwise the column at `objective_column_index`
will be updated.

Note, the weights are negated inside the function so that the objective metabolites 
are seen as reagents/substrates, not products in the reaction equation. 

# Example
```
add_objective!(model, ["met1", "met2"]) # adds a new column with weights = [1,1]
add_objective!(model, ["met1", "met2"]; objective_weights=[0.1, 0.9]) # adds a new column
add_objective!(model, ["met1", "met2"]; objective_weights=[0.1, 0.9], objective_column_index=10) # updates column 10
```
"""
function add_objective!(
    community::CoreModel,
    objective_mets::Vector{String};
    objective_weights = Float64[],
    objective_column_index = 0,
)
    obj_inds = indexin(objective_mets, metabolites(community))
    if isempty(objective_weights)
        objective_weights = repeat([1.0], inner = length(objective_mets))
    end

    if objective_column_index == 0 # needs to be created
        nr, _ = size(community.S)
        objcol = spzeros(nr)
        objcol[obj_inds] .= -objective_weights

        # extend model by one reaction 
        community.S = hcat(community.S, objcol)
        community.xl = [community.xl; 0.0]
        community.xu = [community.xu; 1000.0]
        community.rxns = [community.rxns; "community_biomass"]
        community.c = spzeros(size(community.S, 2))
        community.c[end] = 1.0
    else # only update
        community.S[:, objective_column_index] .= 0.0 # reset
        community.S[obj_inds, objective_column_index] .= -objective_weights
        dropzeros!(community.S)
        community.c = spzeros(size(community.S, 2))
        community.c[objective_column_index] = 1.0
        community.xl[objective_column_index] = 0.0
        community.xu[objective_column_index] = 1000.0
    end
end

add_objective!(model::CoreModel, objective_met::String, objective_weight::Float64) =
    add_objective!(model, [objective_met], [objective_weight])

"""
add_model(
    community::CoreModel,
    model::M,
    exchange_rxn_ids::Vector{String},
    exchange_met_ids::Vector{String};
    species_name="",
    biomass_id=""
    ) where {M<:MetabolicModel}

Add `model` to `community`, which is a pre-existing community model with
`exchange_rxn_ids` and `exchange_met_ids`. The `species_name` is appended to
each reaction and metabolite, see [`join`](@ref). If `biomass_id` is specified
then a biomass metabolite for `model` is also added to the resulting model. The
column corresponding to the `biomass_id` reaction then produces this new biomass
metabolite with unit coefficient. Note, `exchange_rxn_ids` and
`exchange_met_ids` must already exist in the `community` model.

# Example
```
community = add_model(community, model, exchange_rxn_ids, exchange_met_ids; species_name="species_2", biomass_id="BIOMASS_Ecoli_core_w_GAM")
```
"""
function add_model(
    community::CoreModel,
    model::MetabolicModel,
    exchange_rxn_ids::Vector{String},
    exchange_met_ids::Vector{String};
    species_name = "unknown_species",
    biomass_id = nothing,
)::CoreModel

    exchange_met_community_inds = indexin(exchange_met_ids, metabolites(community))
    exchange_rxn_community_inds = indexin(exchange_rxn_ids, reactions(community))
    if any(isnothing.(exchange_met_community_inds)) ||
       any(isnothing.(exchange_rxn_community_inds))
        throw(
            DomainError(
                "exchange metabolite/reaction not found.",
                "Exchange metabolites/reactions must already be contained in the community model.",
            ),
        )
    end

    n_cmodel_rows, n_cmodel_cols = size(stoichiometry(community))
    n_model_rows, n_model_cols = size(stoichiometry(model))
    # A note on the variable names here.Suppose M is some sparse matrix, then I
    # = row indices, J = column indices and V = values at the associated
    # indices. So I[a] = i, J[a]=j and then M[i,j] = V[a]
    Iadd, Jadd, Vadd = findnz(stoichiometry(model))

    # shift to fit into community         
    Iadd .+= n_cmodel_rows
    Jadd .+= n_cmodel_cols

    # when adding a single model not that many reactions, push! okay?
    exchange_rxn_model_inds = indexin(exchange_rxn_ids, reactions(model))
    for i = 1:length(exchange_rxn_ids)
        isnothing(exchange_rxn_model_inds[i]) && continue
        push!(Iadd, exchange_met_community_inds[i]) # already present ex met in community model
        push!(Jadd, n_cmodel_cols + exchange_rxn_model_inds[i]) # new column hence the offset
        push!(Vadd, 1.0)
    end

    biomass_met = 0.0
    if biomass_id != "" # add biomass metabolite
        biomass_rxn = first(indexin([biomass_id], reactions(model)))
        push!(Iadd, n_model_rows + n_cmodel_rows + 1)
        push!(Jadd, biomass_rxn + n_cmodel_cols)
        push!(Vadd, 1.0)
        biomass_met = 1
    end

    n_metabolites_total = n_model_rows + n_cmodel_rows + biomass_met
    n_reactions_total = n_cmodel_cols + n_model_cols

    I, J, V = findnz(stoichiometry(community))
    I = [I; Iadd]
    J = [J; Jadd]
    V = [V; Vadd]
    S = sparse(I, J, V, n_metabolites_total, n_reactions_total)

    # A note on the variables here. The bounds are vectors of upper and lower
    # bounds for each reaction. So lbs = [lb_1, lb_2, lb_i, ...], ubs = [ub_1,
    # ub_2, ub_i, ...] for reaction i. See the bounds function for more
    # information 
    lbsadd, ubsadd = bounds(model)
    lbs, ubs = bounds(community)
    lbs = [lbs; lbsadd]
    ubs = [ubs; ubsadd]

    rxnsadd = "$(species_name)_" .* reactions(model)
    if !isnothing(biomass_id)
        metsadd =
            ["$(species_name)_" .* metabolites(model); "$(species_name)_" * biomass_id]
    else
        metsadd = "$(species_name)_" .* metabolites(model)
    end
    rxns = [reactions(community); rxnsadd]
    mets = [metabolites(community); metsadd]

    # adds to the original community data, could possibly reset?
    I, V = findnz(balance(community))
    b = sparsevec(I, V, n_metabolites_total)
    I, V = findnz(objective(community))
    c = sparsevec(I, V, n_reactions_total)

    return CoreModel(S, b, c, lbs, ubs, rxns, mets)
end

"""
    join_with_exchanges(models::Vector{M}, 
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
m1 = load_model(core_model_path)
m2 = load_model(CoreModel, core_model_path)

boundary_rxn_ids, boundary_met_ids = boundary_reactions_metabolites(m2)
exchange_rxn_ids = filter(startswith("EX_"), boundary_rxn_ids)
exchange_met_ids = filter(endswith("_e"), boundary_met_ids)

biomass_ids = ["BIOMASS_Ecoli_core_w_GAM", "BIOMASS_Ecoli_core_w_GAM"]

community = join_with_exchanges(
    [m1, m2],
    exchange_rxn_ids,
    exchange_met_ids;
    add_biomass_objective = true,
    biomass_ids = biomass_ids,
)
```
"""
function join_with_exchanges(
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

    if length(exchange_met_ids) != length(exchange_rxn_ids)
        throw(
            DomainError(
                "Exchange identifiers are misspecified",
                "The lenghts of the exchange metabolite and reactions are different.",
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
    boundary_reactions_metabolites(model::MetabolicModel)::Tuple{Vector{String}, Vector{String}}

Return `boundary_rxns` and `boundary_mets` given `model`. Boundary reactions
are reactions with a single stoichiometric coefficient in a column in the
stoichiometric matrix, and boundary metabolites are the corresponding row/metabolite
for that column.
"""
function boundary_reactions_metabolites(
    model::MetabolicModel,
)::Tuple{Vector{String},Vector{String}}
    boundary_mets = String[]
    boundary_rxns = String[]
    S = stoichiometry(model)
    rxns = reactions(model)
    mets = metabolites(model)
    for i = 1:size(S, 2)
        n = nnz(S[:, i])
        if n == 1
            j, _ = findnz(S[:, i])
            push!(boundary_mets, mets[first(j)])
            push!(boundary_rxns, rxns[i])
        end
    end
    return boundary_rxns, boundary_mets
end
