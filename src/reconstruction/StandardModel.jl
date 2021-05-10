"""
    add!(model::StandardModel, rxns::Union{Vector{Reaction}, Reaction})

Add `rxn(s)` to `model` if they are not already present based on reaction `id`.
"""
function add!(model::StandardModel, rxns::Vector{Reaction})
    for rxn in rxns
        add!(model, rxn)
    end
end

function add!(model::StandardModel, rxn::Reaction)
    model.reactions[rxn.id] = rxn
end

"""
    add!(model::StandardModel, mets::Union{Vector{Metabolite}, Metabolite})

Add `met(s)` to `model` if they are not already present, based on metabolite `id`.
"""
function add!(model::StandardModel, mets::Vector{Metabolite})
    for met in mets
        add!(model, met)
    end
end

function add!(model::StandardModel, met::Metabolite)
    model.metabolites[met.id] = met
end

"""
    add!(model::StandardModel, genes::Union{Vector{Gene}, Gene})

Add `gene(s)` to `model` if they are not already present based on gene `id`.
"""
function add!(model::StandardModel, genes::Vector{Gene})
    for gene in genes
        add!(model, gene)
    end
end

function add!(model::StandardModel, gene::Gene)
    model.genes[gene.id] = gene
end

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
        push!(all_reactions.args, :(add!($model, r)))
    end
    return all_reactions
end

"""
    rm!(::Type{Reaction}, model::StandardModel, ids::Vector{String})

Remove all reactions with `ids` from `model`.

# Example
rm!(Reaction, model, ["EX_glc__D_e", "fba"])
rm!(Reaction, model, "EX_glc__D_e")
"""
function rm!(::Type{Reaction}, model::StandardModel, ids::Vector{String})
    for id in ids
        rm!(Reaction, model, id)
    end
end

function rm!(::Type{Reaction}, model::StandardModel, id::String)
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
