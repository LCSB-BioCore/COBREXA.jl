# Add and remove reactions

"""
$(TYPEDSIGNATURES)

Add `rxns` to `model` based on reaction `id`.
"""
function add_reactions!(model::ObjectModel, rxns::Vector{Reaction})
    rxn_ids = collect(keys(model.reactions))
    idxs = filter(!isnothing, indexin([r.id for r in rxns], rxn_ids))
    isempty(idxs) ||
        throw(ArgumentError("Duplicated reaction IDs in model: $(rxn_ids[idxs])"))

    for rxn in rxns
        model.reactions[rxn.id] = rxn
    end
end

"""
$(TYPEDSIGNATURES)

Add `rxn` to `model` based on reaction `id`.
"""
add_reaction!(model::ObjectModel, rxn::Reaction) = add_reactions!(model, [rxn])

"""
$(TYPEDSIGNATURES)

Add `rxns` to `model` and return a shallow copied version of the model.
"""
function add_reactions(model::ObjectModel, rxns::Vector{Reaction})
    m = copy(model)

    m.reactions = copy(m.reactions)
    for rxn in rxns
        m.reactions[rxn.id] = rxn
    end

    m
end

"""
$(TYPEDSIGNATURES)

Add `rxn` to `model`, and return a shallow copied version of the model.
"""
add_reaction(model::ObjectModel, rxn::Reaction) = add_reactions(model, [rxn])

@_remove_fn reaction ObjectModel String inplace begin
    if !(reaction_id in variables(model))
        @models_log @info "Reaction $reaction_id not found in model."
    else
        delete!(model.reactions, reaction_id)
    end
    nothing
end

@_remove_fn reaction ObjectModel String inplace plural begin
    remove_reaction!.(Ref(model), reaction_ids)
    nothing
end

@_remove_fn reaction ObjectModel String begin
    remove_reactions(model, [reaction_id])
end

@_remove_fn reaction ObjectModel String plural begin
    n = copy(model)
    n.reactions = copy(model.reactions)
    remove_reactions!(n, reaction_ids)
    return n
end

# Add and remove metabolites

"""
$(TYPEDSIGNATURES)

Add `mets` to `model` based on metabolite `id`.
"""
function add_metabolites!(model::ObjectModel, mets::Vector{Metabolite})
    met_ids = collect(keys(model.metabolites))
    idxs = filter(!isnothing, indexin([m.id for m in mets], met_ids))
    isempty(idxs) ||
        throw(ArgumentError("Duplicated metabolite IDs in model: $(met_ids[idxs])"))

    for met in mets
        model.metabolites[met.id] = met
    end
end

"""
$(TYPEDSIGNATURES)

Add `met` to `model` based on metabolite `id`.
"""
add_metabolite!(model::ObjectModel, met::Metabolite) = add_metabolites!(model, [met])

"""
$(TYPEDSIGNATURES)

Add `mets` to `model` and return a shallow copied version of the model.
"""
function add_metabolites(model::ObjectModel, mets::Vector{Metabolite})
    m = copy(model)

    m.metabolites = copy(m.metabolites)
    for met in mets
        m.metabolites[met.id] = met
    end

    m
end

"""
$(TYPEDSIGNATURES)

Add `met` to `model` and return a shallow copied version of the model.
"""
add_metabolite(model::ObjectModel, met::Metabolite) = add_metabolites(model, [met])

@_remove_fn metabolite ObjectModel String inplace begin
    remove_metabolites!(model, [metabolite_id])
end

@_remove_fn metabolite ObjectModel String inplace plural begin
    !all(in.(metabolite_ids, Ref(metabolites(model)))) &&
        @models_log @info "Some metabolites not found in model."
    remove_reactions!(
        model,
        [
            rid for (rid, rn) in model.reactions if
            any(haskey.(Ref(rn.metabolites), metabolite_ids))
        ],
    )
    delete!.(Ref(model.metabolites), metabolite_ids)
    nothing
end

@_remove_fn metabolite ObjectModel String begin
    remove_metabolites(model, [metabolite_id])
end

@_remove_fn metabolite ObjectModel String plural begin
    n = copy(model)
    n.reactions = copy(model.reactions)
    n.metabolites = copy(model.metabolites)
    remove_metabolites!(n, metabolite_ids)
    return n
end

# Add and remove genes

"""
$(TYPEDSIGNATURES)

Add `genes` to `model` based on gene `id`.
"""
function add_genes!(model::ObjectModel, genes::Vector{Gene})
    gene_ids = collect(keys(model.genes))
    idxs = filter(!isnothing, indexin([g.id for g in genes], gene_ids))
    isempty(idxs) || throw(ArgumentError("Duplicated gene IDs in model: $(gene_ids[idxs])"))

    for gene in genes
        model.genes[gene.id] = gene
    end
end

"""
$(TYPEDSIGNATURES)

Add `gene` to `model` based on gene `id`.
"""
add_gene!(model::ObjectModel, gene::Gene) = add_genes!(model, [gene])

"""
$(TYPEDSIGNATURES)

Add `gns` to `model` and return a shallow copied version of the model.
"""
function add_genes(model::ObjectModel, genes::Vector{Gene})
    m = copy(model)

    m.genes = copy(m.genes)
    for gn in genes
        m.genes[gn.id] = gn
    end

    m
end

"""
$(TYPEDSIGNATURES)

Add `gene` to `model` and return a shallow copied version of the model.
"""
add_gene(model::ObjectModel, gene::Gene) = add_genes(model, [gene])

"""
$(TYPEDSIGNATURES)

Remove all genes with `ids` from `model`. If `knockout_reactions` is true, then also
constrain reactions that require the genes to function to carry zero flux.

# Example
```
remove_genes!(model, ["g1", "g2"])
```
"""
function remove_genes!(
    model::ObjectModel,
    gids::Vector{String};
    knockout_reactions::Bool = false,
)
    if knockout_reactions
        rm_reactions = String[]
        for (rid, r) in model.reactions
            if !isnothing(r.gene_associations) && all(
                any(in.(gids, Ref(conjunction))) for
                conjunction in reaction_gene_associations(model, rid)
            )
                push!(rm_reactions, rid)
            end
        end
        pop!.(Ref(model.reactions), rm_reactions)
    end
    pop!.(Ref(model.genes), gids)
    nothing
end

"""
$(TYPEDSIGNATURES)

Remove gene with `id` from `model`. If `knockout_reactions` is true, then also
constrain reactions that require the genes to function to carry zero flux.

# Example
```
remove_gene!(model, "g1")
```
"""
remove_gene!(model::ObjectModel, gid::String; knockout_reactions::Bool = false) =
    remove_genes!(model, [gid]; knockout_reactions = knockout_reactions)

# Change reaction bounds

@_change_bounds_fn ObjectModel String inplace begin
    isnothing(lower_bound) || (model.reactions[rxn_id].lower_bound = lower_bound)
    isnothing(upper_bound) || (model.reactions[rxn_id].upper_bound = upper_bound)
    nothing
end

@_change_bounds_fn ObjectModel String inplace plural begin
    for (i, l, u) in zip(rxn_ids, lower_bounds, upper_bounds)
        change_bound!(model, i, lower_bound = l, upper_bound = u)
    end
end

@_change_bounds_fn ObjectModel String begin
    change_bounds(
        model,
        [rxn_id],
        lower_bounds = [lower_bound],
        upper_bounds = [upper_bound],
    )
end

@_change_bounds_fn ObjectModel String plural begin
    n = copy(model)
    n.reactions = copy(model.reactions)
    for rid in rxn_ids
        n.reactions[rid] = copy(model.reactions[rid])
        for field in fieldnames(typeof(model.reactions[rid]))
            setfield!(n.reactions[rid], field, getfield(model.reactions[rid], field))
        end
    end
    change_bounds!(n, rxn_ids; lower_bounds, upper_bounds)
    return n
end

# Change gene product bounds

"""
$(TYPEDSIGNATURES)

Changes the `product_lower_bound` or `product_upper_bound` for the
[`Gene`][(ref)s listed in `gids` in the `model`, in place.
"""
function change_gene_product_bounds!(
    model::ObjectModel,
    gids::Vector{String};
    lower_bounds = fill(nothing, length(gids)),
    upper_bounds = fill(nothing, length(gids)),
)
    for (i, gid) in enumerate(gids)

        isnothing(lower_bounds[i]) || begin
            (model.genes[gid].product_lower_bound = lower_bounds[i])
        end

        isnothing(upper_bounds[i]) || begin
            (model.genes[gid].product_upper_bound = upper_bounds[i])
        end
    end
end

"""
$(TYPEDSIGNATURES)

Changes the `product_lower_bound` or `product_upper_bound` for the
[`Gene`][(ref) `gid` in the `model`, in place.
"""
function change_gene_product_bound!(
    model::ObjectModel,
    gid::String;
    lower_bound = nothing,
    upper_bound = nothing,
)
    change_gene_product_bounds!(
        model,
        [gid];
        lower_bounds = [lower_bound],
        upper_bounds = [upper_bound],
    )
end

"""
$(TYPEDSIGNATURES)

Variant of [`change_gene_product_bounds!`](@ref) that does not modify the
original model, but makes a shallow copy with the modification included.
"""
function change_gene_product_bounds(
    model::ObjectModel,
    gids::Vector{String};
    lower_bounds = fill(nothing, length(gids)),
    upper_bounds = fill(nothing, length(gids)),
)
    m = copy(model)
    m.genes = copy(model.genes)
    for gid in gids
        m.genes[gid] = copy(model.genes[gid])
        for field in fieldnames(typeof(m.genes[gid]))
            setfield!(m.genes[gid], field, getfield(model.genes[gid], field))
        end
    end
    change_gene_product_bounds!(m, gids; lower_bounds, upper_bounds)
    m
end

"""
$(TYPEDSIGNATURES)

Variant of [`change_gene_product_bound!`](@ref) that does not modify the
original model, but makes a shallow copy with the modification included.
"""
function change_gene_product_bound(
    model::ObjectModel,
    gid::String;
    lower_bound = nothing,
    upper_bound = nothing,
)
    change_gene_product_bounds(
        model,
        [gid];
        lower_bounds = [lower_bound],
        upper_bounds = [upper_bound],
    )
end

# Change objective 

"""
$(TYPEDSIGNATURES)

Change the objective for `model` to reaction(s) with `rxn_ids`, optionally
specifying their `weights`. By default, assume equal weights. If no objective
exists in model, sets objective.
"""
function change_objective!(
    model::ObjectModel,
    rxn_ids::Vector{String};
    weights = ones(length(rxn_ids)),
)
    all(!haskey(model.reactions, rid) for rid in rxn_ids) &&
        throw(DomainError(rxn_ids, "Some reaction ids were not found in model."))
    model.objective = Dict(rxn_ids .=> weights)
    nothing
end

"""
$(TYPEDSIGNATURES)

Variant of [`change_objective!`](@ref) that sets a single `rxn_id` as the
objective weight with `weight` (defaults to 1.0).
"""
change_objective!(model::ObjectModel, rxn_id::String; weight::Float64 = 1.0) =
    change_objective!(model, [rxn_id]; weights = [weight])

"""
$(TYPEDSIGNATURES)

Variant of [`change_objective!`](@ref) that does not modify the original model,
but makes a shallow copy with the modification included.
"""
function change_objective(
    model::ObjectModel,
    rxn_ids::Vector{String};
    weights = ones(length(rxn_ids)),
)
    m = copy(model)
    m.objective = copy(model.objective)
    change_objective!(m, rxn_ids; weights)
    m
end

"""
$(TYPEDSIGNATURES)

Variant of [`change_objective!`](@ref) that does not modify the original model,
but makes a shallow copy with the modification included.
"""
function change_objective(model::ObjectModel, rxn_id::String; weight::Float64 = 1.0)
    m = copy(model)
    m.objective = copy(model.objective)
    change_objective!(m, rxn_id; weight)
    m
end

# Add and remove biomass metabolite

"""
$(TYPEDSIGNATURES)

Add a biomass metabolite called `biomass_metabolite_id` with stoichiometry 1 to
the biomass reaction, called `biomass_rxn_id` in `model`. Changes the model in
place. Does not check if the model already has a biomass metabolite.
"""
function add_biomass_metabolite!(
    model::ObjectModel,
    biomass_rxn_id::String;
    biomass_metabolite_id = "biomass",
)
    model.reactions[biomass_rxn_id].metabolites[biomass_metabolite_id] = 1.0
    add_metabolite!(model, Metabolite(biomass_metabolite_id))
end

"""
$(TYPEDSIGNATURES)

Variant of [`add_biomass_metabolite!`](@ref) that does not modify the original
model, but makes a shallow copy with the modification included. Does not check
if the model already has a biomass metabolite.
"""
function add_biomass_metabolite(
    model::ObjectModel,
    biomass_rxn_id::String;
    biomass_metabolite_id = "biomass",
)
    m = copy(model)
    m.metabolites = copy(m.metabolites)
    m.metabolites[biomass_metabolite_id] = Metabolite(biomass_metabolite_id)

    m.reactions = copy(model.reactions)
    m.reactions[biomass_rxn_id] = copy(model.reactions[biomass_rxn_id])
    m.reactions[biomass_rxn_id].metabolites =
        copy(model.reactions[biomass_rxn_id].metabolites)
    m.reactions[biomass_rxn_id].metabolites[biomass_metabolite_id] = 1.0

    m
end

"""
$(TYPEDSIGNATURES)

Remove a biomass metabolite called `biomass_metabolite_id` from
the biomass reaction, called `biomass_rxn_id` in `model`. 
"""
function remove_biomass_metabolite!(
    model::ObjectModel,
    biomass_rxn_id::String;
    biomass_metabolite_id = "biomass",
)
    delete!(model.reactions[biomass_rxn_id].metabolites, biomass_metabolite_id)
    remove_metabolite!(model, biomass_metabolite_id)
end

"""
$(TYPEDSIGNATURES)

Variant of [`remove_biomass_metabolite!`](@ref) that does not modify the original
model, but makes a shallow copy with the modification included.
"""
function remove_biomass_metabolite(
    model::ObjectModel,
    biomass_rxn_id::String;
    biomass_metabolite_id = "biomass",
)
    m = copy(model)
    m.metabolites = copy(m.metabolites)
    remove_metabolite!(m, biomass_metabolite_id)

    m.reactions = copy(model.reactions)
    m.reactions[biomass_rxn_id] = copy(model.reactions[biomass_rxn_id])
    m.reactions[biomass_rxn_id].metabolites =
        copy(model.reactions[biomass_rxn_id].metabolites)
    delete!(m.reactions[biomass_rxn_id].metabolites, biomass_metabolite_id)

    m
end

# Add virtual ribosome for constrained allocation applications

"""
$(TYPEDSIGNATURES)

To `biomass_rxn_id` in `model`, add a pseudo-isozyme that approximates the
effect ribosome synthesis has on growth.

# Bacterial growth law models
Numerous experimental studies have shown that during steady state growth the
cellular density of E. coli (and probably other bacteria) is constant. As a
consequence of this, growth law models typically assume that the total proteome
capacity (mass fraction of protein in the cell) is limited. Further, it has been
shown experimentally that the ribosomal protein content of a bacterial cell
increases with faster growth rate. Ribosomes are used to make proteins (and
themselves), leading to a trade-off: faster growth requires more ribosomes to
make more enzymes to grow, but this reduces the amount of proteome space "left
over" for biosynthetic enzymes. See Mori, Matteo, et al. "Constrained allocation
flux balance analysis." PLoS computational biology 12.6 (2016) for more details.

# Implementation
By adding a protein cost to biomass synthesis this effect can be simulated. The
parameter `weight` needs to be estimated from data, but acts like a turnover
number. Lower `weight` means more ribosome is required for growth (`ribosome =
growth/weight`). The molar mass of the ribosome is `1`.

# Note
1. This modifications makes the underlying biomass reaction unidirectional.
2. The `virtualribosome_id` defaults to `virtualribosome` and must be manually included
   in any capacity bound later used in enzyme constrained models.
3. This modification also adds a pseudogene called `virtualribosome_id`.

The pseudo-isozyme acts like a regular gene product,
```
ribosome = weight * biomass_flux
```
"""
function add_virtualribosome!(
    model::ObjectModel,
    biomass_rxn_id::String,
    weight::Float64;
    virtualribosome_id = "virtualribosome",
)
    # ensure unidirectional
    model.reactions[biomass_rxn_id].lower_bound = 0.0
    model.reactions[biomass_rxn_id].upper_bound = constants.default_reaction_bound

    # add ribosome kinetics
    model.reactions[biomass_rxn_id].gene_associations = [
        Isozyme(
            kcat_forward = weight,
            kcat_backward = 0.0,
            gene_product_stoichiometry = Dict(virtualribosome_id => 1.0),
        ),
    ]

    # add ribosome gene
    model.genes[virtualribosome_id] =
        Gene(id = virtualribosome_id, product_molar_mass = 1.0)

    nothing
end

"""
$(TYPEDSIGNATURES)

Variant of [`add_virtualribosome!`](@ref) that does not modify the original
model, but makes a shallow copy with the modification included.
"""
function add_virtualribosome(
    model::ObjectModel,
    biomass_rxn_id::String,
    weight::Float64;
    virtualribosome_id = "virtualribosome",
)
    m = copy(model)
    m.reactions = copy(model.reactions)
    m.reactions[biomass_rxn_id] = copy(model.reactions[biomass_rxn_id])
    m.genes = copy(model.genes)

    add_virtualribosome!(m, biomass_rxn_id, weight; virtualribosome_id)

    m
end

# Add, change, remove isozymes

"""
$(TYPEDSIGNATURES)

Add `isozymes` to `rxn_id` in `model`.
"""
add_isozymes!(model::ObjectModel, rxn_id::String, isozymes::Vector{Isozyme}) =
    add_isozymes!(model, [rxn_id], [isozymes])

"""
$(TYPEDSIGNATURES)

For each reaction in `rxn_ids`, add the corresponding isozymes in
`isozymes_vector` to `model`.
"""
function add_isozymes!(
    model::ObjectModel,
    rxn_ids::Vector{String},
    isozymes_vector::Vector{Vector{Isozyme}},
)
    rxn_ids = keys(model.reactions)
    idxs = filter(!isnothing, indexin(rxns, rxn_ids))
    isempty(idxs) ||
        throw(ArgumentError("Duplicated reaction IDs in model: $(rxn_ids[idxs])"))
    isnothing(model.reactions[rxn_id].gene_associations) ||
        throw(ArgumentError("$rxn_id already has isozymes."))

    for (rid, isozymes) in zip(rxn_id_vector, isozymes_vector)
        add_isozymes!(model, rid, isozymes)
    end
end

"""
$(TYPEDSIGNATURES)

Variant of [`add_isozymes!`](@ref) that returns a copied model instead of
modifying the input.
"""
function add_isozymes(model::ObjectModel, rxn_id::String, isozymes::Vector{Isozyme})

    m = copy(model)
    m.reactions = copy(model.reactions)

    m.reactions[rxn_id] = copy(model.reactions[rxn_id])
    if !isnothing(model.reactions[rxn_id].gene_associations)
        m.reactions[rxn_id].gene_associations =
            copy(model.reactions[rxn_id].gene_associations)
    end
    m.reactions[rxn_id].gene_associations = isozymes

    m
end

"""
$(TYPEDSIGNATURES)

For each pair of `isozymes` and `rxn_id` in `isozymes_vector` and
`rxn_id_vector`, call [`add_isozymes`](@ref) to add the isozymes to the
`model`.
"""
function add_isozymes(
    model::ObjectModel,
    rxn_id_vector::Vector{String},
    isozymes_vector::Vector{Vector{Isozyme}},
)

    m = copy(model)
    m.reactions = copy(model.reactions)

    for (rxn_id, isozymes) in zip(rxn_id_vector, isozymes_vector)

        m.reactions[rxn_id] = copy(model.reactions[rxn_id])
        if !isnothing(model.reactions[rxn_id].gene_associations)
            m.reactions[rxn_id].gene_associations =
                copy(model.reactions[rxn_id].gene_associations)
        end
        m.reactions[rxn_id].gene_associations = isozymes

    end

    m
end
