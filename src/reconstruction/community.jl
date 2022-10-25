"""
$(TYPEDSIGNATURES)

Add an objective column to the `community` model with optional id `objective_id`. Supply a
dictionary mapping the string names of the objective metabolites to their weights in
`objective_mets_weights`. Note, the weights are negated inside the function so that positive
weights are seen as reagents/substrates, NOT products in the reaction equation.

# Example
```
add_community_objective!(model, Dict("met1"=>1.0, "met2"=>2.0))
```

See also: [`update_community_objective!`](@ref)
"""
function add_community_objective!(
    community::CoreModel,
    objective_mets_weights::Dict{String,Float64};
    objective_id = "community_biomass",
)
    obj_inds = indexin(keys(objective_mets_weights), metabolites(community))
    nothing in obj_inds && throw(ArgumentError("Some objective metabolite(s) not found."))

    nr, _ = size(community.S)
    objcol = spzeros(nr)
    objcol[obj_inds] .= -collect(values(objective_mets_weights))

    # extend model by one reaction
    community.S = hcat(community.S, objcol)
    community.xl = [community.xl; 0.0]
    community.xu = [community.xu; constants.default_reaction_bound]
    community.rxns = [community.rxns; objective_id]
    community.c = spzeros(size(community.S, 2))
    community.c[end] = 1.0

    return nothing # stop annoying return value
end

"""
$(TYPEDSIGNATURES)

Variant of [`add_community_objective!`] that takes a `StandardModel` community model as input.
"""
function add_community_objective!(
    community::StandardModel,
    objective_mets_weights::Dict{String,Float64};
    objective_id = "community_biomass",
)
    nothing in indexin(keys(objective_mets_weights), metabolites(community)) &&
        throw(ArgumentError("Some objective metabolite(s) not found."))

    rdict = Dict(k => -float(v) for (k, v) in objective_mets_weights)
    rxn = Reaction(objective_id)
    rxn.metabolites = rdict
    rxn.lower_bound = 0.0
    rxn.upper_bound = constants.default_reaction_bound
    rxn.objective_coefficient = 1.0
    community.reactions[rxn.id] = rxn

    return nothing # stop annoying return value
end

"""
$(TYPEDSIGNATURES)

Update the weights for the objective column with id `objective_id` in `community` using
`objective_mets_weights`, which maps metabolite ids to weights. The current weights are
reset to 0 before being updated to the supplied weights. Note, the weights are negated
inside the function so that the objective metabolites are seen as reagents/substrates, NOT
products in the reaction equation.

# Example
```
update_community_objective!(model, "community_biomass", Dict("met1"=>1.0, "met2"=>2.0))
```

See also: [`add_community_objective!`](@ref)
"""
function update_community_objective!(
    community::CoreModel,
    objective_id::String,
    objective_mets_weights::Dict{String,Float64},
)
    obj_inds = indexin(keys(objective_mets_weights), metabolites(community))
    nothing in obj_inds && throw(ArgumentError("Some objective metabolite(s) not found."))

    objective_column_index = first(indexin([objective_id], reactions(community)))
    community.S[:, objective_column_index] .= 0.0 # reset
    community.S[obj_inds, objective_column_index] .=
        -collect(values(objective_mets_weights))
    dropzeros!(community.S)
    community.c = spzeros(size(community.S, 2))
    community.c[objective_column_index] = 1.0
    community.xl[objective_column_index] = 0.0
    community.xu[objective_column_index] = constants.default_reaction_bound

    return nothing # stop annoying return value
end

"""
$(TYPEDSIGNATURES)

Variant of [`update_community_objective!`] that takes a `StandardModel` community model as input.
"""
function update_community_objective!(
    community::StandardModel,
    objective_id::String,
    objective_mets_weights::Dict{String,Float64},
)
    delete!(community.reactions, objective_id)
    add_community_objective!(community, objective_mets_weights; objective_id = objective_id)

    return nothing # stop annoying return value
end

"""
$(TYPEDSIGNATURES)

Return a `CoreModel` representing the community model of `models` joined through their
exchange reactions and metabolites in the dictionary `exchange_rxn_mets`, which maps
exchange reactions to their associated metabolite. These exchange reactions and metabolites
link model metabolites to environmental metabolites and reactions. Optionally specify
`model_names` to append a specific name to each reaction and metabolite of an organism for
easier reference (default is `species_i` for each model index i in `models`). Note, the
bounds of the environmental variables are all set to zero. Thus, to run a simulation you
need to constrain them appropriately. All the other bounds are inherited from the models
used to construct the community model.

If `biomass_ids` is supplied, then a community model is returned that has an extra reaction
added to the end of the stoichiometric matrix (last column) that can be assigned as the
objective reaction. It also creates biomass "metabolites" that can be used in this objective
reaction. In the returned mode, these biomass metabolites are produced by the reaction
corresponding to `biomass_ids` in each model respectively. Note, this reaction is
unspecified, further action needs to be taken to specify it, e.g. assign weights to the last
column of the stoichiometric matrix in the rows corresponding to the biomass metabolites.

To further clarify how this `join` works. Suppose you have 2 organisms with stoichiometric
matrices S₁ and S₂ and you want to link them with `exchange_rxn_mets = Dict(er₁ => em₁, er₂
=> em₂, er₃ => em₃, ...)`. Then a new community stoichiometric matrix is constructed that
looks like this:
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
The exchange reactions in each model get linked to environmental metabolites, `emᵢ`, and
these get linked to environmental exchanges, `erᵢ`. These `erᵢ` behave like normal single
organism exchange reactions. When `biomass_ids` are supplied, each model's biomass reaction
produces a pseudo-metabolite (`bmᵢ`). These can be weighted in column `b`, called the
`community_biomass` reaction in the community model, if desired. Refer to the tutorial if
this is unclear.

# Example
```
m1 = load_model(core_model_path)
m2 = load_model(CoreModel, core_model_path)

# need to list ALL the exchanges that will form part of the entire model
exchange_rxn_mets = Dict(k => first(keys(reaction_stoichiometry(m1, ex_rxn)))
    for filter(looks_like_exchange_reaction, reactions(m1)))

biomass_ids = ["BIOMASS_Ecoli_core_w_GAM", "BIOMASS_Ecoli_core_w_GAM"]

community = join_with_exchanges(
    CoreModel,
    [m1, m2],
    exchange_rxn_mets;
    biomass_ids = biomass_ids,
)
```
"""
function join_with_exchanges(
    ::Type{CoreModel},
    models::Vector{M},
    exchange_rxn_mets::Dict{String,String};
    biomass_ids = String[],
    model_names = String[],
) where {M<:AbstractMetabolicModel}

    exchange_rxn_ids = keys(exchange_rxn_mets)
    exchange_met_ids = values(exchange_rxn_mets)
    add_biomass_objective = isempty(biomass_ids) ? false : true

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
        length(exchange_met_ids) +
        nnz_add
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

    _reaction_offsets = cumsum(reaction_lengths)
    _metabolite_offsets = cumsum(metabolite_lengths)
    for i = 1:length(models)
        species = isempty(model_names) ? "species_$(i)_" : model_names[i] * "_"
        tlbs, tubs = bounds(models[i])
        lbs[reaction_offset[i]+1:_reaction_offsets[i]] .= tlbs
        ubs[reaction_offset[i]+1:_reaction_offsets[i]] .= tubs
        rxns[reaction_offset[i]+1:_reaction_offsets[i]] = species .* reactions(models[i])
        mets[metabolite_offset[i]+1:_metabolite_offsets[i]] =
            species .* metabolites(models[i])
    end
    mets[_metabolite_offsets[end]+1:_metabolite_offsets[end]+length(exchange_met_ids)] .=
        exchange_met_ids
    rxns[_reaction_offsets[end]+1:_reaction_offsets[end]+length(exchange_rxn_ids)] .=
        exchange_rxn_ids

    if add_biomass_objective
        rxns[end] = "community_biomass"
        for i = 1:length(models)
            species = isempty(model_names) ? "species_$(i)_" : model_names[i] * "_"
            mets[end-length(biomass_ids)+i] = species .* biomass_ids[i]
        end
    end

    return CoreModel(S, spzeros(size(S, 1)), spzeros(size(S, 2)), lbs, ubs, rxns, mets)
end

"""
$(TYPEDSIGNATURES)

A variant of [`join_with_exchanges`](@ref) that returns a `StandardModel`.
"""
function join_with_exchanges(
    ::Type{StandardModel},
    models::Vector{M},
    exchange_rxn_mets::Dict{String,String};
    biomass_ids = [],
    model_names = [],
)::StandardModel where {M<:AbstractMetabolicModel}

    community = StandardModel()
    rxns = OrderedDict{String,Reaction}()
    mets = OrderedDict{String,Metabolite}()
    genes = OrderedDict{String,Gene}()
    sizehint!(rxns, sum(n_reactions(m) for m in models) + length(keys(exchange_rxn_mets)))
    sizehint!(
        mets,
        sum(n_metabolites(m) for m in models) + length(values(exchange_rxn_mets)),
    )
    sizehint!(genes, sum(n_genes(m) for m in models))

    for (i, model) in enumerate(models)
        species = isempty(model_names) ? "species_$(i)" : model_names[i] # underscore gets added in add_model_with_exchanges!
        biomass_id = isempty(biomass_ids) ? nothing : biomass_ids[i]
        add_model_with_exchanges!(
            community,
            model,
            exchange_rxn_mets;
            model_name = species,
            biomass_id = biomass_id,
        )
    end

    if !isempty(biomass_ids)
        bm = Dict{String,Float64}() # unassigned
        community.reactions["community_biomass"] =
            Reaction("community_biomass", bm, :forward)
    end

    # Add environmental exchange reactions and metabolites.TODO: add annotation details
    for (rid, mid) in exchange_rxn_mets
        community.reactions[rid] = Reaction(rid)
        community.reactions[rid].metabolites = Dict{String,Float64}(mid => -1.0)
        community.reactions[rid].lower_bound = 0.0
        community.reactions[rid].upper_bound = 0.0
        community.metabolites[mid] = Metabolite(mid)
        community.metabolites[mid].id = mid
    end

    return community
end

"""
$(TYPEDSIGNATURES)

Add `model` to `community`, which is a pre-existing community model with exchange reactions
and metabolites in the dictionary `exchange_rxn_mets`. The `model_name` is appended to each
reaction and metabolite, see [`join_with_exchanges`](@ref). If `biomass_id` is specified
then a biomass metabolite for `model` is also added to the resulting model. The column
corresponding to the `biomass_id` reaction then produces this new biomass metabolite with
unit coefficient. The exchange reactions and metabolites in `exchange_rxn_mets` must already
exist in `community`. Always returns a new community model because it is more efficient than
resizing all the matrices.

No in-place variant for `CoreModel`s exists yet.

# Example
```
community = add_model_with_exchanges(community,
    model,
    exchange_rxn_mets;
    model_name="species_2",
    biomass_id="BIOMASS_Ecoli_core_w_GAM")
```
"""
function add_model_with_exchanges(
    community::CoreModel,
    model::AbstractMetabolicModel,
    exchange_rxn_mets::Dict{String,String};
    model_name = "unknown_species",
    biomass_id = nothing,
)

    exchange_rxn_ids = keys(exchange_rxn_mets)
    exchange_met_ids = values(exchange_rxn_mets)
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
    exchange_rxn_model_idxs = indexin(exchange_rxn_ids, reactions(model))
    for i = 1:length(exchange_rxn_ids)
        isnothing(exchange_rxn_model_idxs[i]) && continue
        push!(Iadd, exchange_met_community_inds[i]) # already present ex met in community model
        push!(Jadd, n_cmodel_cols + exchange_rxn_model_idxs[i]) # new column hence the offset
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

    rxnsadd = "$(model_name)_" .* reactions(model)
    if !isnothing(biomass_id)
        metsadd = ["$(model_name)_" .* metabolites(model); "$(model_name)_" * biomass_id]
    else
        metsadd = "$(model_name)_" .* metabolites(model)
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
$(TYPEDSIGNATURES)

The `StandardModel` variant of [`add_model_with_exchanges`](@ref), but is in-place.
"""
function add_model_with_exchanges!(
    community::StandardModel,
    model::AbstractMetabolicModel,
    exchange_rxn_mets::Dict{String,String};
    model_name = "unknown_species",
    biomass_id = nothing,
)
    stdm = model isa StandardModel ? deepcopy(model) : convert(StandardModel, model)
    model_name = model_name * "_"

    for met in values(stdm.metabolites)
        met.id = model_name * met.id
        community.metabolites[met.id] = met
    end

    for rxn in values(stdm.reactions)
        # rename reaction string
        rxn.metabolites = Dict(model_name * k => v for (k, v) in rxn.metabolites)
        # change direction of exchange
        if rxn.id in keys(exchange_rxn_mets)
            rxn.metabolites[model_name*exchange_rxn_mets[rxn.id]] = 1.0 # exchanges should be negative originally. TODO: test if they are?
        end
        # add biomass metabolite if applicable
        if rxn.id == biomass_id
            rxn.metabolites[model_name*rxn.id] = 1.0 # produces one biomass
            community.metabolites[model_name*rxn.id] = Metabolite(model_name * rxn.id)
        end
        # add environmental connection
        if rxn.id in keys(exchange_rxn_mets)
            rxn.metabolites[exchange_rxn_mets[rxn.id]] = -1.0
        end
        # reset objective
        rxn.objective_coefficient = 0.0
        rxn.id = model_name * rxn.id
        # add to community model
        community.reactions[rxn.id] = rxn
    end

    for gene in values(stdm.genes)
        gene.id = model_name * gene.id
        community.genes[gene.id] = gene
    end

    return nothing
end

"""
$(TYPEDSIGNATURES)

The `StandardModel` variant of [`add_model_with_exchanges`](@ref). Makes a deepcopy of
`community` and calls the inplace variant of this function on that copy.
"""
function add_model_with_exchanges(
    community::StandardModel,
    model::AbstractMetabolicModel,
    exchange_rxn_mets::Dict{String,String};
    model_name = "unknown_species",
    biomass_id = nothing,
)
    new_comm = deepcopy(community)
    add_model_with_exchanges!(
        new_comm,
        model,
        exchange_rxn_mets;
        model_name = model_name,
        biomass_id = biomass_id,
    )
    return new_comm
end
