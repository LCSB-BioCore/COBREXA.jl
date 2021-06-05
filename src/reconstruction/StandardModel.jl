"""
    add_reactions!(model::StandardModel, rxns::Vector{Reaction})

Add `rxns` to `model` based on reaction `id`.
"""
function add_reactions!(model::StandardModel, rxns::Vector{Reaction})
    for rxn in rxns
        model.reactions[rxn.id] = rxn
    end
end

add_reactions!(model::StandardModel, rxn::Reaction) = add_reactions!(model, [rxn])

"""
    add_metabolites!(model::StandardModel, mets::Vector{Metabolite})

Add `mets` to `model` based on metabolite `id`.
"""
function add_metabolites!(model::StandardModel, mets::Vector{Metabolite})
    for met in mets
        model.metabolites[met.id] = met
    end
end

add_metabolites!(model::StandardModel, met::Metabolite) =  add_metabolite!(model, [met])
    
"""
    add_genes!(model::StandardModel, genes::Vector{Gene})

Add `genes` to `model` based on gene `id`.
"""
function add_genes!(model::StandardModel, genes::Vector{Gene})
    for gene in genes
        model.genes[gene.id] = gene
    end
end

add_genes!(model::StandardModel, gene::Gene) = add_genes!(model, [gene]) 

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
        push!(all_reactions.args, :(add_reactions!($model, r)))
    end
    return all_reactions
end

"""
    remove_reactions!(model::StandardModel, ids::Vector{String})

Remove all reactions with `ids` from `model`.

# Example
remove_reactions!(model, ["EX_glc__D_e", "fba"])
"""
function remove_reactions!(model::StandardModel, ids::Vector{String})
    for id in ids
        rm!(Reaction, model, id)
    end
end

"""
    remove_reaction!(model::StandardModel, id::String)

Remove reaction with `id` from `model`.

# Example
remove_reaction!(model, "EX_glc__D_e")
"""
function remove_reaction!(model::StandardModel, id::String)
    rxn = pop!(model.reactions, id)
end

"""
    rm!(::Type{Metabolite}, model::StandardModel, ids::Vector{String})

Remove all metabolites with `ids` from `model`.
Warning, this could leave the model inconsistent, e.g. a reaction might
require the deleted metabolite.

# Example
rm!(Metabolite, model, ["atp_c", "adp_c"])
rm!(Metabolite, model, "atp_c")
"""
function rm!(::Type{Metabolite}, model::StandardModel, ids::Vector{String})
    for id in ids
        rm!(Metabolite, model, id)
    end
end

function rm!(::Type{Metabolite}, model::StandardModel, id::String)
    delete!(model.metabolites, id)
end

"""
    rm!(::Type{Gene}, model::StandardModel, ids::Vector{String})

Remove all genes with `ids` from `model`.

# Example
rm!(Gene, model, ["g1", "g2"])
rm!(Gene, model, "g1")
"""
rm!(::Type{Gene}, model::StandardModel, gid::String; knockout_reactions::Bool = false) =
    rm!(Gene, model, [gid]; knockout_reactions = knockout_reactions)

function rm!(
    ::Type{Gene},
    model::StandardModel,
    gids::Vector{String};
    knockout_reactions::Bool = false,
)
    if knockout_reactions
        rm_reactions = String[]
        for (rid, r) in model.reactions
            if !isnothing(r.grr) &&
               all([any(in.(gids, Ref(conjunction))) for conjunction in r.grr])
                push!(rm_reactions, rid)
            end
        end
        delete!.(Ref(model.reactions), rm_reactions)
    end
    delete!.(Ref(model.genes), gids)
end

function set_bound(model::StandardModel, reaction_id::String; ub, lb)
    reaction = model.reactions[reaction_id]
    reaction.lb = lb
    reaction.ub = ub
end
