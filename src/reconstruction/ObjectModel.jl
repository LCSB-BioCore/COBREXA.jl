# Add and remove reactions

"""
$(TYPEDSIGNATURES)

Plural variant of [`add_reaction!`](@ref).
"""
function add_reactions!(model::ObjectModel, rxns::Vector{Reaction})
    throw_argerror_if_key_found(model, :reactions, rxns)
    for rxn in rxns
        model.reactions[rxn.id] = rxn
    end
end

"""
$(TYPEDSIGNATURES)

Add `rxn` to `model` based on reaction `id` if the `id` is not already in the
model.
"""
add_reaction!(model::ObjectModel, rxn::Reaction) = add_reactions!(model, [rxn])

"""
$(TYPEDSIGNATURES)

Plural variant of [`add_reaction`](@ref).
"""
function add_reactions(model::ObjectModel, rxns::Vector{Reaction})
    m = copy(model)
    m.reactions = copy(m.reactions)
    add_reactions!(m, rxns)
    m
end

"""
$(TYPEDSIGNATURES)

Return a shallow copied version of `model` with `rxn` added if it's ID is not
already present in the model.
"""
add_reaction(model::ObjectModel, rxn::Reaction) = add_reactions(model, [rxn])

@_remove_fn reaction ObjectModel String inplace plural begin
    throw_argerror_if_key_missing(model, :reactions, reaction_ids)
    delete!.(Ref(model.reactions), reaction_ids)
    nothing
end

@_remove_fn reaction ObjectModel String inplace begin
    remove_reactions!(model, [reaction_id])
    nothing
end

@_remove_fn reaction ObjectModel String plural begin
    n = copy(model)
    n.reactions = copy(model.reactions)
    remove_reactions!(n, reaction_ids)
    return n
end

@_remove_fn reaction ObjectModel String begin
    remove_reactions(model, [reaction_id])
end

# Add and remove metabolites

"""
$(TYPEDSIGNATURES)

Plural variant of [`add_metabolite!`](@ref).
"""
function add_metabolites!(model::ObjectModel, mets::Vector{Metabolite})
    throw_argerror_if_key_found(model, :metabolites, mets)
    for met in mets
        model.metabolites[met.id] = met
    end
end

"""
$(TYPEDSIGNATURES)

Add `met` to `model` based on metabolite `id` if the `id` is not already in the
model.
"""
add_metabolite!(model::ObjectModel, met::Metabolite) = add_metabolites!(model, [met])

"""
$(TYPEDSIGNATURES)

Plural variant of [`add_metabolite`](@ref).
"""
function add_metabolites(model::ObjectModel, mets::Vector{Metabolite})
    m = copy(model)
    m.metabolites = copy(m.metabolites)
    add_metabolites!(m, mets)
    m
end

"""
$(TYPEDSIGNATURES)

Return a shallow copied version of the `model` with `met` added.  Only adds
`met` if its ID is not already present in the model.
"""
add_metabolite(model::ObjectModel, met::Metabolite) = add_metabolites(model, [met])

@_remove_fn metabolite ObjectModel String inplace plural begin
    throw_argerror_if_key_missing(model, :metabolites, metabolite_ids)
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

@_remove_fn metabolite ObjectModel String inplace begin
    remove_metabolites!(model, [metabolite_id])
end

@_remove_fn metabolite ObjectModel String plural begin
    n = copy(model)
    n.reactions = copy(model.reactions)
    n.metabolites = copy(model.metabolites)
    remove_metabolites!(n, metabolite_ids)
    return n
end

@_remove_fn metabolite ObjectModel String begin
    remove_metabolites(model, [metabolite_id])
end

# Add and remove genes

"""
$(TYPEDSIGNATURES)

Plural variant of [`add_gene!`](@ref).
"""
function add_genes!(model::ObjectModel, genes::Vector{Gene})
    throw_argerror_if_key_found(model, :genes, genes)
    for gene in genes
        model.genes[gene.id] = gene
    end
end

"""
$(TYPEDSIGNATURES)

Add `gene` to `model` based on gene `id` if the `id` is not already in the
model.
"""
add_gene!(model::ObjectModel, gene::Gene) = add_genes!(model, [gene])

"""
$(TYPEDSIGNATURES)

Plural variant of [`add_gene`](@ref).
"""
function add_genes(model::ObjectModel, genes::Vector{Gene})
    m = copy(model)
    m.genes = copy(m.genes)
    add_genes!(m, genes)
    m
end

"""
$(TYPEDSIGNATURES)

Return a shallow copied version of the `model` with added `gene`. Only adds the
`gene` if its ID is not already present in the model.
"""
add_gene(model::ObjectModel, gene::Gene) = add_genes(model, [gene])

"""
$(TYPEDSIGNATURES)

Plural variant of [`remove_gene!`](@ref).
"""
function remove_genes!(
    model::ObjectModel,
    gids::Vector{String};
    knockout_reactions::Bool = false,
)
    throw_argerror_if_key_missing(model, :genes, gids)
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
        delete!.(Ref(model.reactions), rm_reactions)
    end
    delete!.(Ref(model.genes), gids)
    nothing
end

"""
$(TYPEDSIGNATURES)

Remove gene with `id` from `model`. If `knockout_reactions` is true, then also
constrain reactions that require the genes to function to carry zero flux.
"""
remove_gene!(model::ObjectModel, gid::String; knockout_reactions::Bool = false) =
    remove_genes!(model, [gid]; knockout_reactions = knockout_reactions)

# Change reaction bounds

@_change_bounds_fn ObjectModel String inplace begin
    change_bounds!(
        model,
        [rxn_id];
        lower_bounds = [lower_bound],
        upper_bounds = [upper_bound],
    )
end

@_change_bounds_fn ObjectModel String inplace plural begin
    throw_argerror_if_key_missing(model, :reactions, rxn_ids)
    for (rxn_id, lower, upper) in zip(rxn_ids, lower_bounds, upper_bounds)
        isnothing(lower) || (model.reactions[rxn_id].lower_bound = lower)
        isnothing(upper) || (model.reactions[rxn_id].upper_bound = upper)
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
    throw_argerror_if_key_missing(model, :reactions, rxn_ids)
    m = copy(model)
    m.reactions = copy(model.reactions)
    for (rid, lower, upper) in zip(rxn_ids, lower_bounds, upper_bounds)
        m.reactions[rid] = copy(model.reactions[rid])
        for field in fieldnames(typeof(model.reactions[rid]))
            setfield!(m.reactions[rid], field, getfield(model.reactions[rid], field))
        end
        isnothing(lower) || (m.reactions[rxn_id].lower_bound = lower)
        isnothing(upper) || (m.reactions[rxn_id].upper_bound = upper)
    end
    return m
end

# Change gene product bounds

"""
$(TYPEDSIGNATURES)

Plural variant of [`change_gene_product_bound!`](@ref).
"""
function change_gene_product_bounds!(
    model::ObjectModel,
    gids::Vector{String};
    lower_bounds = fill(nothing, length(gids)),
    upper_bounds = fill(nothing, length(gids)),
)
    throw_argerror_if_key_missing(model, :genes, gids)
    for (gid, lower, upper) in zip(gids, lower_bounds, upper_bounds)
        isnothing(lower) || (model.genes[gid].product_lower_bound = lower)
        isnothing(upper) || (model.genes[gid].product_upper_bound = upper)
    end
end

"""
$(TYPEDSIGNATURES)

Changes the `product_lower_bound` or `product_upper_bound` for the
[`Gene`][(ref) `gid` in the `model`, in place. If either `lower_bound` or
`upper_bound` is `nothing`, then that bound is not changed.
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

Plural variant of [`change_gene_product_bound`](@ref).
"""
function change_gene_product_bounds(
    model::ObjectModel,
    gids::Vector{String};
    lower_bounds = fill(nothing, length(gids)),
    upper_bounds = fill(nothing, length(gids)),
)
    throw_argerror_if_key_missing(model, :genes, gids)
    m = copy(model)
    m.genes = copy(model.genes)
    for (gid, lower, upper) in zip(gids, lower_bounds, upper_bounds)
        m.genes[gid] = copy(model.genes[gid])
        for field in fieldnames(typeof(model.genes[gid]))
            setfield!(m.genes[gid], field, getfield(model.genes[gid], field))
        end
        isnothing(lower) || (model.genes[gid].product_lower_bound = lower)
        isnothing(upper) || (model.genes[gid].product_upper_bound = upper)
    end
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
    haskey(model.reactions, biomass_rxn_id) ||
        throw(ArgumentError("$biomass_rxn_id not found in model."))
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
    haskey(model.reactions, biomass_rxn_id) ||
        throw(ArgumentError("$biomass_rxn_id not found in model."))

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
    haskey(model.reactions, biomass_rxn_id) ||
        throw(ArgumentError("$biomass_rxn_id not found in model."))
    haskey(model.reaction[biomass_rxn_id], biomass_metabolite_id) ||
        throw(ArgumentError("$biomass_metabolite_id not found in $biomass_rxn_id."))

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
    haskey(model.reactions, biomass_rxn_id) ||
        throw(ArgumentError("$biomass_rxn_id not found in model."))
    haskey(model.reaction[biomass_rxn_id], biomass_metabolite_id) ||
        throw(ArgumentError("$biomass_metabolite_id not found in $biomass_rxn_id."))

    m = copy(model)
    m.metabolites = copy(m.metabolites)
    delete!(m.metabolites, biomass_metabolite_id)

    m.reactions = copy(model.reactions)
    m.reactions[biomass_rxn_id] = copy(model.reactions[biomass_rxn_id])
    m.reactions[biomass_rxn_id].metabolites =
        copy(model.reactions[biomass_rxn_id].metabolites)
    delete!(m.reactions[biomass_rxn_id].metabolites, biomass_metabolite_id)

    m
end

# Add virtual ribosome (assume no model already has one)

"""
$(TYPEDSIGNATURES)

To `biomass_rxn_id` in `model`, add a pseudo-isozyme and associated gene that
approximates the effect ribosome synthesis has on growth.

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
This modification makes the underlying biomass reaction unidirectional. The
`virtualribosome_id` defaults to `virtualribosome`, and corresponds to a
pseudogene called `virtualribosome_id`. The parameter `weight` needs to be
estimated from data, but acts like a turnover number. Lower `weight` means more
ribosome is required for growth (`ribosome = growth/weight`). The molar mass of
the ribosome is `1`. The pseudo-isozyme acts like a regular gene product,
```
ribosome = weight * biomass_flux
```
when simulating enzyme constrained models.
"""
function add_virtualribosome!(
    model::ObjectModel,
    biomass_rxn_id::String,
    weight::Float64;
    virtualribosome_id = "virtualribosome",
)
    haskey(model.reactions, biomass_rxn_id) ||
        throw(ArgumentError("$biomass_rxn_id not found in model."))
    isnothing(model.reactions[biomass_rxn_id]) ||
        throw(ArgumentError("$biomass_rxn_id already has isozymes associated to it."))
    haskey(model.genes, virtualribosome_id) ||
        throw(ArgumentError("$virtualribosome_id already found in model."))

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
    haskey(model.reactions, biomass_rxn_id) ||
        throw(ArgumentError("$biomass_rxn_id not found in model."))
    isnothing(model.reactions[biomass_rxn_id]) ||
        throw(ArgumentError("$biomass_rxn_id already has isozymes associated to it."))
    haskey(model.genes, virtualribosome_id) ||
        throw(ArgumentError("$virtualribosome_id already found in model."))

    m = copy(model)
    m.reactions = copy(model.reactions)
    m.reactions[biomass_rxn_id] = copy(model.reactions[biomass_rxn_id])
    m.genes = copy(model.genes)

    add_virtualribosome!(m, biomass_rxn_id, weight; virtualribosome_id)

    m
end

# Add, remove isozymes (no change because the order isozymes may appear in is not constant across models)

"""
$(TYPEDSIGNATURES)

Plural variant of [`add_isozymes!`](@ref).
"""
function add_isozymes!(
    model::ObjectModel,
    rids::Vector{String},
    isozymes_vector::Vector{Vector{Isozyme}},
)
    throw_argerror_if_key_missing(model, :reactions, rids)
    throw_argerror_if_isozymes_found(model, rids)

    for (rid, isozymes) in zip(rids, isozymes_vector)
        model.reactions[rid].gene_associations = isozymes
    end
end

"""
$(TYPEDSIGNATURES)

Add `isozymes` to `rid` in `model`. Only allowed if `rid` does not have isozymes assigned to it.
"""
add_isozymes!(model::ObjectModel, rid::String, isozymes::Vector{Isozyme}) =
    add_isozymes!(model, [rid], [isozymes])

"""
$(TYPEDSIGNATURES)

Plural variant of [`add_isozymes`](@ref).
"""
function add_isozymes(
    model::ObjectModel,
    rids::Vector{String},
    isozymes_vector::Vector{Vector{Isozyme}},
)
    throw_argerror_if_key_missing(model, :reactions, rids)
    throw_argerror_if_isozymes_found(model, rids)

    m = copy(model)
    m.reactions = copy(model.reactions)

    for (rid, isozymes) in zip(rids, isozymes_vector)
        m.reactions[rid] = copy(model.reactions[rid])
        if !isnothing(model.reactions[rid].gene_associations)
            m.reactions[rid].gene_associations =
                copy(model.reactions[rid].gene_associations)
        end
        m.reactions[rid].gene_associations = isozymes
    end
    m
end

"""
$(TYPEDSIGNATURES)

Variant of [`add_isozymes!`](@ref) that returns a copied model instead of
modifying the input.
"""
add_isozymes(model::ObjectModel, rid::String, isozymes::Vector{Isozyme}) =
    add_isozymes(model, [rid], [isozymes])

"""
$(TYPEDSIGNATURES)

Plural variant of [`remove_isozyme!`](@ref).
"""
function remove_isozymes!(model::ObjectModel, rids::Vector{String})
    throw_argerror_if_key_missing(model, :reactions, rids)
    for rid in rids
        model.reactions[rid].gene_associations = nothing
    end
end

"""
$(TYPEDSIGNATURES)

Remove all isozymes from `rid` in `model`.
"""
remove_isozyme!(model::ObjectModel, rid::String) = remove_isozymes!(model, [rid])

"""
$(TYPEDSIGNATURES)

Plural variant of [`remove_isozyme`](@ref).
"""
function remove_isozymes!(model::ObjectModel, rids::Vector{String})
    throw_argerror_if_key_missing(model, :reactions, rids)

    m = copy(model)
    m.reactions = copy(model.reactions)
    for rid in rids
        m.reactions[rid] = copy(model.reactions[rid])
        for field in fieldnames(typeof(model.reactions[rid]))
            setfield!(m.reactions[rid], field, getfield(model.reactions[rid], field))
        end
        model.reactions[rid].gene_associations = nothing
    end
end

"""
$(TYPEDSIGNATURES)

Return a shallow copy of `model` with all isozymes of reaction `rid` removed.
"""
remove_isozyme(model::ObjectModel, rid::String) = remove_isozymes(model, [rid])
