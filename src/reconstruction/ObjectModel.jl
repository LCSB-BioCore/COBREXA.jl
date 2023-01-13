"""
$(TYPEDSIGNATURES)

Add `rxns` to `model` based on reaction `id`.
"""
function add_reactions!(model::ObjectModel, rxns::Vector{Reaction})
    for rxn in rxns
        model.reactions[rxn.id] = rxn
    end
    nothing
end

"""
$(TYPEDSIGNATURES)

Add `rxn` to `model` based on reaction `id`.
"""
add_reaction!(model::ObjectModel, rxn::Reaction) = add_reactions!(model, [rxn])

"""
$(TYPEDSIGNATURES)

Add `mets` to `model` based on metabolite `id`.
"""
function add_metabolites!(model::ObjectModel, mets::Vector{Metabolite})
    for met in mets
        model.metabolites[met.id] = met
    end
    nothing
end

"""
$(TYPEDSIGNATURES)

Add `met` to `model` based on metabolite `id`.
"""
add_metabolite!(model::ObjectModel, met::Metabolite) = add_metabolites!(model, [met])

"""
$(TYPEDSIGNATURES)

Add `genes` to `model` based on gene `id`.
"""
function add_genes!(model::ObjectModel, genes::Vector{Gene})
    for gene in genes
        model.genes[gene.id] = gene
    end
    nothing
end

"""
$(TYPEDSIGNATURES)

Add `gene` to `model` based on gene `id`.
"""
add_gene!(model::ObjectModel, gene::Gene) = add_genes!(model, [gene])

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


@_change_bounds_fn ObjectModel String inplace begin
    isnothing(lower) || (model.reactions[rxn_id].lower_bound = lower)
    isnothing(upper) || (model.reactions[rxn_id].upper_bound = upper)
    nothing
end

@_change_bounds_fn ObjectModel String inplace plural begin
    for (i, l, u) in zip(rxn_ids, lower, upper)
        change_bound!(model, i, lower = l, upper = u)
    end
end

@_change_bounds_fn ObjectModel String begin
    change_bounds(model, [rxn_id], lower = [lower], upper = [upper])
end

@_change_bounds_fn ObjectModel String plural begin
    n = copy(model)
    n.reactions = copy(model.reactions)
    for i in rxn_ids
        n.reactions[i] = copy(n.reactions[i])
    end
    change_bounds!(n, rxn_ids, lower = lower, upper = upper)
    return n
end

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

"""
$(TYPEDSIGNATURES)

Change the objective for `model` to reaction(s) with `rxn_ids`, optionally specifying their `weights`. By default,
assume equal weights. If no objective exists in model, sets objective.
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

change_objective!(model::ObjectModel, rxn_id::String) = change_objective!(model, [rxn_id])

"""
$(TYPEDSIGNATURES)

Add a biomass metabolite called `biomass_metabolite_id` with stoichiometry 1 to
the biomass reaction, called `biomass_rxn_id` in `model`. Changes the model in
place.
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
model, but makes a shallow copy with the modification included.
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
2. The `pseudoribosome_id` defaults to `pseudoribosome` and must be manually included
   in any capacity bound later used in enzyme constrained models.
3. This modification also adds a pseudogene called `pseudoribosome_id`.

The pseudo-isozyme acts like a regular gene product,
```
ribosome = weight * biomass_flux
```
"""
function add_pseudoribosome!(
    model::ObjectModel,
    biomass_rxn_id::String,
    weight::Float64;
    pseudoribosome_id = "pseudoribosome",
)
    # ensure unidirectional
    model.reactions[biomass_rxn_id].lower_bound = 0.0
    model.reactions[biomass_rxn_id].upper_bound = constants.default_reaction_bound

    # add ribosome kinetics
    model.reactions[biomass_rxn_id].gene_associations = [
        Isozyme(
            kcat_forward = weight,
            kcat_backward = 0.0,
            gene_product_stoichiometry = Dict(pseudoribosome_id => 1.0),
        ),
    ]

    # add ribosome gene
    model.genes[pseudoribosome_id] = Gene(id = pseudoribosome_id, product_molar_mass = 1.0)

    nothing
end

"""
$(TYPEDSIGNATURES)

Variant of [`add_pseudoribosome!`](@ref) that does not modify the original
model, but makes a shallow copy with the modification included.
"""
function add_pseudoribosome(
    model::ObjectModel,
    biomass_rxn_id::String,
    weight::Float64;
    pseudoribosome_id = "pseudoribosome",
)
    m = copy(model)

    # unidirectional biomass
    m.reactions = copy(model.reactions)
    m.reactions[biomass_rxn_id] = copy(model.reactions[biomass_rxn_id])
    m.reactions[biomass_rxn_id].lower_bound = 0.0
    m.reactions[biomass_rxn_id].upper_bound = constants.default_reaction_bound

    # add ribosome kinetics
    if !isnothing(model.reactions[biomass_rxn_id].gene_associations)
        m.reactions[biomass_rxn_id].gene_associations =
            copy(model.reactions[biomass_rxn_id].gene_associations)
    end
    m.reactions[biomass_rxn_id].gene_associations = [
        Isozyme(
            kcat_forward = weight,
            kcat_backward = 0.0,
            gene_product_stoichiometry = Dict(pseudoribosome_id => 1.0),
        ),
    ]

    # add ribosome gene
    m.genes = copy(model.genes)
    m.genes[pseudoribosome_id] = Gene(id = pseudoribosome_id, product_molar_mass = 1.0)

    m
end

"""
$(TYPEDSIGNATURES)

Add `isozymes` to `rxn_id` in `model`. Overwrites the currently stored isozymes. Assumes genes
are already in `model`.
"""
add_isozymes!(
    model::ObjectModel,
    rxn_id::String,
    isozymes::Vector{Isozyme},
) = model.reactions[rxn_id].gene_associations = isozymes

"""
$(TYPEDSIGNATURES)

For each pair of `isozymes` and `rxn_id` in `isozymes_vector` and
`rxn_id_vector`, call [`add_isozymes!`](@ref) to add the isozymes to the
`model`.
"""
function add_isozymes!(
    model::ObjectModel,
    rxn_id_vector::Vector{String},
    isozymes_vector::Vector{Vector{Isozyme}},
)
    for (rid, isozymes) in zip(rxn_id_vector, isozymes_vector)
        add_isozymes!(model, rid, isozymes)
    end
end

"""
$(TYPEDSIGNATURES)

Variant of [`add_isozymes!`](@ref) that returns a copied model instead of
modifying the input. 
"""
function add_isozymes(
    model::ObjectModel,
    rxn_id::String,
    isozymes::Vector{Isozyme},
)

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
