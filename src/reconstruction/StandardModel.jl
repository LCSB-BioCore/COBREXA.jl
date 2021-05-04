function _map_reaction_to_genes!(model::StandardModel, rxn::Reaction)
    if !isnothing(rxn.grr)
        for gene_array in rxn.grr
            for gene in gene_array
                push!(model.genes[gene].associated_reactions, rxn.id)
            end
        end
    end
end

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
    _map_reaction_to_genes!(model, rxn)
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
    if !isnothing(rxn.grr)
        for gene_array in rxn.grr
            for gene_id in gene_array
                delete!(model.genes[gene_id].associated_reactions, id)
            end
        end
    end
end

"""
    rm!(::Type{Metabolite}, model::StandardModel, ids::Vector{String})

Remove all metabolites with `ids` from `model`.

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
function rm!(
    ::Type{Gene},
    model::StandardModel,
    ids::Vector{String};
    knockout::Bool = false,
)
    for id in ids
        rm!(Gene, model, id, knockout = knockout)
    end
end

function rm!(::Type{Gene}, model::StandardModel, id::String; knockout::Bool = false)
    if knockout
        # Dev note: the three nested for loops are inefficiency. However:
        # - gene_ids (user input) will be probably only very few items
        # - model.genes[gene_id].reactions are just a few reactions (most genes don't code for a lot of reactions)
        # - reaction.grr also should only hold few items (reactions aren't coded by many different combinations of genes)
        # Let's avoid premature optimization for now and see if anyone ever has problems with this
        for reaction_id in gene_associated_reactions(model, id)
            reaction = get(model.reactions, reaction_id, nothing)
            if isnothing(reaction)
                return nothing
            end
            # AND inside the gene_array, so destroy as soon as one is missing
            reaction.grr = filter(x -> !any(occursin.(id, x)), reaction.grr)

            # OR outside, so all have to be deleted for the reaction to be deleted
            if length(reaction.grr) == 0
                rm!(Reaction, model, reaction.id)
            end
        end
    end
    delete!(model.genes, id)
end

function set_bound(model::StandardModel, reaction_id::String; ub, lb)
    reaction = model.reactions[reaction_id]
    reaction.lb = lb
    reaction.ub = ub
end
