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
    if model[rxn] == -1
        push!(model.reactions, rxn)
    else
        @warn "$(rxn.id) already present in model."
    end
    return nothing
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
    if model[met] == -1
        push!(model.metabolites, met)
    else
        @warn "$(met.id) already present in model."
    end
    return nothing
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
    if model[gene] == -1
        push!(model.genes, gene)
    else
        @warn "$(gene.id) already present in model."
    end
    return nothing
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
function rm!(model::StandardModel, rxns::Union{Vector{Reaction},Reaction})
    new_rxn_list = Reaction[]
    for r in model.reactions
        if typeof(rxns) == Reaction
            if rxns.id != r.id
                push!(new_rxn_list, r)
            end
        else
            if !(r.id in [rr.id for rr in rxns])
                push!(new_rxn_list, r)
            end
        end
    end
    model.reactions = new_rxn_list
    return nothing
end

"""
    rm!(model::StandardModel, mets::Union{Vector{Metabolite}, Metabolite})

Remove `met(s)` from `model` based on metabolite `id`.
"""
function rm!(model::StandardModel, mets::Union{Vector{Metabolite},Metabolite})
    new_met_list = Metabolite[]
    for m in model.metabolites
        if typeof(mets) == Metabolite
            if mets.id != m.id
                push!(new_met_list, m)
            end
        else
            if !(m.id in [mm.id for mm in mets])
                push!(new_met_list, m)
            end
        end
    end
    model.metabolites = new_met_list
    return nothing
end

"""
    rm!(model::StandardModel, genes::Union{Vector{Gene}, Gene})

Remove `gene(s)` from `model` based on gene `id`.
"""
function rm!(model::StandardModel, genes::Union{Vector{Gene},Gene})
    new_gene_list = Gene[]
    for g in model.genes
        if typeof(genes) == Gene
            if genes.id != g.id
                push!(new_gene_list, g)
            end
        else
            if !(g.id in [gg.id for gg in genes])
                push!(new_gene_list, g)
            end
        end
    end
    model.genes = new_gene_list
    return nothing
end

"""
    fix_model!(model::StandardModel)

Inspect metabolites and genes of `model` relative to the reactions of `model`.
Remove genes or metabolites that are not used in the reactions.
Add genes or metabolites that are not present in `model.genes` or `model.metabolites` but are used in `model.reactions`.
Everything is based on the `id` of metabolites and genes, thus is it possible for a metabolite or gene to be duplicated but using
a different `id`. Other functions are provided to help identify these cases.

See also: [`check_duplicate_annotations`](@ref), [`check_duplicate_reaction`](@ref), [`check_same_formula`](@ref).
"""
function fix_model!(model::StandardModel)
    rxn_mets = Metabolite[] # list of metabolites used in reactions
    for rxn in model.reactions
        for met in keys(rxn.metabolites)
            if rxn_mets[met] == -1
                push!(rxn_mets, met)
            end
        end
    end

    rxn_genes = Gene[] # list of genes used in reactions
    for rxn in model.reactions
        for gene_list in rxn.grr # for [] in [[]]
            for gene in gene_list
                if rxn_genes[gene] == -1
                    push!(rxn_genes, gene)
                end
            end
        end
    end

    model_gene_ids = [x.id for x in model.genes]
    model_mets_ids = [x.id for x in model.metabolites]
    rxn_gene_ids = [x.id for x in rxn_genes]
    rxn_mets_ids = [x.id for x in rxn_mets]

    extra_genes = setdiff(model_gene_ids, rxn_gene_ids)
    !isempty(extra_genes) && rm!(model, [findfirst(model.genes, x) for x in extra_genes])
    extra_mets = setdiff(model_mets_ids, rxn_mets_ids)
    !isempty(extra_mets) &&
        rm!(model, [findfirst(model.metabolites, x) for x in extra_mets])

    missing_genes = setdiff(rxn_gene_ids, model_gene_ids)
    for mg in missing_genes
        g = findfirst(rxn_genes, mg)
        add!(model, g)
    end
    missing_mets = setdiff(rxn_mets_ids, model_mets_ids)
    for mm in missing_mets
        m = findfirst(rxn_mets, mm)
        add!(model, m)
    end
end
