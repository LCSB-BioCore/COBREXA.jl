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
    v1, ∅ ⟶ A, 0, 500
    v2, A ⟷ B, -500
    v3, B ⟶ ∅
end
```
"""
macro add_reactions!(model::Symbol, ex::Expr)
    model = esc(model)
    all_reactions = Expr(:block)
    for line in MacroTools.striplines(ex).args
        args = line.args
        id = esc(args[1])
        name = String(args[1])
        reaction = esc(args[2])
        push!(all_reactions.args, :($id = $reaction))
        push!(all_reactions.args, :($id.id = $name))
        if length(args) == 3
            lb = args[3]
            push!(all_reactions.args, :($id.lb = $lb))
        elseif length(args) == 4
            lb = args[3]
            ub = args[4]
            push!(all_reactions.args, :($id.lb = $lb))
            push!(all_reactions.args, :($id.ub = $ub))
        end
        push!(all_reactions.args, :(add!($model, $id)))
    end
    return all_reactions
end

"""
    rm!(model::StandardModel, rxns::Union{Vector{Reaction}, Reaction})

Remove all `rxn(s)` from `model` if the `id`s match those in `rxns`.
"""
function rm!(model::StandardModel, rxns::Vector{Reaction})
    for rxn in rxns
        rm!(model, rxn)
    end
end

function rm!(model::StandardModel, rxn::Reaction)
    delete!(model.reactions, rxn.id)
end

"""
    rm!(model::StandardModel, mets::Union{Vector{Metabolite}, Metabolite})

Remove `met(s)` from `model` based on metabolite `id`.
"""
function rm!(model::StandardModel, mets::Vector{Metabolite})
    for m in mets
        rm!(model, m)
    end
end

function rm!(model::StandardModel, met::Metabolite)
    delete!(model.metabolites, met.id)
end

"""
    rm!(model::StandardModel, genes::Union{Vector{Gene}, Gene})

Remove `gene(s)` from `model` based on gene `id`.
"""
function rm!(model::StandardModel, genes::Vector{Gene})
    for gene in genes
        rm!(model, gene)
    end
end

function rm!(model::StandardModel, gene::Gene)
    delete!(model.genes, gene.id)
end
