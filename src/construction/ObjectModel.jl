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

Shortcut to add multiple reactions and their lower and upper bounds

Call variants
-------------
```
@add_reactions! model begin
    reaction_name, reaction
end

@add_reactions! model begin
    reaction_name, reaction, lower_bound
end

@add_reactions! model begin
    reaction_name, reaction, lower_bound, upper_bound
end
```

Examples
--------
```
@add_reactions! model begin
    "v1", nothing → A, 0, 500
    "v2", A ↔ B + C, -500
    "v3", B + C → nothing
end
```
"""
macro add_reactions!(model::Symbol, ex::Expr)
    model = esc(model)
    all_reactions = Expr(:block)
    for line in MacroTools.striplines(ex).args
        args = line.args
        id = esc(args[1])
        reaction = esc(args[2])
        push!(all_reactions.args, :(r = $reaction))
        push!(all_reactions.args, :(r.id = $id))
        if length(args) == 3
            lb = args[3]
            push!(all_reactions.args, :(r.lower_bound = $lb))
        elseif length(args) == 4
            lb = args[3]
            ub = args[4]
            push!(all_reactions.args, :(r.lower_bound = $lb))
            push!(all_reactions.args, :(r.upper_bound = $ub))
        end
        push!(all_reactions.args, :(add_reaction!($model, r)))
    end
    return all_reactions
end

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
                conjunction in reaction_gene_association(model, rid)
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
    if !(reaction_id in reactions(model))
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
