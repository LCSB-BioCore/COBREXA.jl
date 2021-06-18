"""
    add_reactions!(model::StandardModel, rxns::Vector{Reaction})

Add `rxns` to `model` based on reaction `id`.
"""
function add_reactions!(model::StandardModel, rxns::Vector{Reaction})
    for rxn in rxns
        model.reactions[rxn.id] = rxn
    end
end

"""
    add_reaction!(model::StandardModel, rxn::Reaction)

Add `rxn` to `model` based on reaction `id`.
"""
add_reaction!(model::StandardModel, rxn::Reaction) = add_reactions!(model, [rxn])

"""
    add_metabolites!(model::StandardModel, mets::Vector{Metabolite})

Add `mets` to `model` based on metabolite `id`.
"""
function add_metabolites!(model::StandardModel, mets::Vector{Metabolite})
    for met in mets
        model.metabolites[met.id] = met
    end
end

"""
    add_metabolite!(model::StandardModel, met::Metabolite)

Add `met` to `model` based on metabolite `id`.
"""
add_metabolite!(model::StandardModel, met::Metabolite) = add_metabolites!(model, [met])

"""
    add_genes!(model::StandardModel, genes::Vector{Gene})

Add `genes` to `model` based on gene `id`.
"""
function add_genes!(model::StandardModel, genes::Vector{Gene})
    for gene in genes
        model.genes[gene.id] = gene
    end
end

"""
    add_gene!(model::StandardModel, genes::Gene)

Add `gene` to `model` based on gene `id`.
"""
add_gene!(model::StandardModel, gene::Gene) = add_genes!(model, [gene])

"""
    @add_reactions!(model::Symbol, ex::Expr)

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
    "v1", nothing ⟶ A, 0, 500
    "v2", A ⟷ B + C, -500
    "v3", B + C ⟶ nothing
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
            push!(all_reactions.args, :(r.lb = $lb))
        elseif length(args) == 4
            lb = args[3]
            ub = args[4]
            push!(all_reactions.args, :(r.lb = $lb))
            push!(all_reactions.args, :(r.ub = $ub))
        end
        push!(all_reactions.args, :(add_reaction!($model, r)))
    end
    return all_reactions
end

"""
    remove_reactions!(model::StandardModel, ids::Vector{String})

Remove all reactions with `ids` from `model`. Note, may result in orphan metabolites.

# Example
```
remove_reactions!(model, ["EX_glc__D_e", "fba"])
```
"""
function remove_reactions!(model::StandardModel, ids::Vector{String})
    pop!.(Ref(model.reactions), ids)
end

"""
    remove_reaction!(model::StandardModel, id::String)

Remove reaction with `id` from `model`. Note, may result in orphan metabolites.

# Example
```
remove_reaction!(model, "EX_glc__D_e")
```
"""
remove_reaction!(model::StandardModel, id::String) = remove_reactions!(model, [id])

"""
    remove_metabolites!(model::StandardModel, ids::Vector{String})

Remove all metabolites with `ids` from `model`.
Warning, this could leave the model inconsistent, e.g. a reaction might
require the deleted metabolite, in which case analysis functions will error.

# Example
```
remove_metabolites!(model, ["atp_c", "adp_c"])
```
"""
function remove_metabolites!(model::StandardModel, ids::Vector{String})
    pop!.(Ref(model.metabolites), ids)
end

"""
    remove_metabolite!(model::StandardModel, id::String)

Remove metabolite with `id` from `model`.
Warning, this could leave the model inconsistent, e.g. a reaction might
require the deleted metabolite, in which case analysis functions will error.

# Example
```
remove_metabolite!(model, "atp_c")
```
"""
remove_metabolite!(model::StandardModel, id::String) = remove_metabolites!(model, [id])

"""
    remove_genes!(
        model::StandardModel,
        ids::Vector{String};
        knockout_reactions::Bool = false,
    )

Remove all genes with `ids` from `model`. If `knockout_reactions` is true, then also 
constrain reactions that require the genes to function to carry zero flux.

# Example
```
remove_genes!(model, ["g1", "g2"])
```
"""
function remove_genes!(
    model::StandardModel,
    gids::Vector{String};
    knockout_reactions::Bool = false,
)
    if knockout_reactions
        rm_reactions = String[]
        for (rid, r) in model.reactions
            if !isnothing(r.grr) &&
               all(any(in.(gids, Ref(conjunction))) for conjunction in r.grr)
                push!(rm_reactions, rid)
            end
        end
        pop!.(Ref(model.reactions), rm_reactions)
    end
    pop!.(Ref(model.genes), gids)
end

"""
    remove_gene!(
        model::StandardModel,
        id::Vector{String};
        knockout_reactions::Bool = false,
    )

Remove gene with `id` from `model`. If `knockout_reactions` is true, then also 
constrain reactions that require the genes to function to carry zero flux.

# Example
```
remove_gene!(model, "g1")
```
"""
remove_gene!(model::StandardModel, gid::String; knockout_reactions::Bool = false) =
    remove_genes!(model, [gid]; knockout_reactions = knockout_reactions)

@doc @_change_bound_s_bang("StandardModel", "rxn_id", "String", "\"PFL\"", false, true)
function change_bound!(
    model::StandardModel,
    reaction_id::String;
    lower_bound = -_constants.default_reaction_bound,
    upper_bound = _constants.default_reaction_bound,
)
    reaction = model.reactions[reaction_id]
    reaction.lb = float(lower_bound)
    reaction.ub = float(upper_bound)
    return nothing # so that nothing gets printed
end

@doc @_change_bound_s_bang("StandardModel", "rxn_ids", "Vector{String}", "[\"PFL\", \"FBA\"]", true, true)
function change_bounds!(
    model::StandardModel,
    reaction_ids::Vector{String};
    lower_bounds = fill(-_constants.default_reaction_bound, length(reaction_ids)),
    upper_bounds = fill(_constants.default_reaction_bound, length(reaction_ids)),
)
    for (rid, lb, ub) in zip(reaction_ids, lower_bounds, upper_bounds)
        change_bound!(model, rid; lower_bound = lb, upper_bound = ub)
    end
end

@doc @_change_bound_s_bang("StandardModel", "rxn_id", "String", "\"PFL\"", false, false)
function change_bound(
    model::StandardModel,
    reaction_id::String;
    lower_bound = -_constants.default_reaction_bound,
    upper_bound = _constants.default_reaction_bound,
)
    m = copy(model)
    m.reactions = copy(model.reactions)
    r = m.reactions[reaction_id] = copy(model.reactions[reaction_id])
    r.lb = float(lower_bound)
    r.ub = float(upper_bound)
    return m
end

@doc @_change_bound_s_bang("StandardModel", "rxn_ids", "Vector{String}", "[\"PFL\", \"FBA\"]", true, false)
function change_bounds(
    model::StandardModel,
    reaction_ids::Vector{String};
    lower_bounds = fill(-_constants.default_reaction_bound, length(reaction_ids)),
    upper_bounds = fill(constants.default_reaction_bound, length(reaction_ids)),
)
    m = copy(model)
    m.reactions = copy(model.reactions)
    for (rid, lb, ub) in zip(reaction_ids, lower_bounds, upper_bounds)
        r = m.reactions[rid] = copy(model.reactions[rid])
        r.lb = float(lb)
        r.ub = float(ub)
    end
    return m
end
