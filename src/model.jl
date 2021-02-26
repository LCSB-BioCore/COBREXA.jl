"""
Model struct of a constraint based metabolic model (mutable).

# Fields
````
id :: String
reactions :: Array{Reaction, 1}
metabolites :: Array{Metabolite, 1}
genes :: Array{Gene, 1}
````
"""
mutable struct Model
    id :: String
    reactions :: Array{Reaction, 1} 
    metabolites :: Array{Metabolite, 1}
    genes :: Array{Gene, 1}
end

"""
Model()

Empty model constructor.
"""
function Model()
    CobraTools.Model("blank", Array{Reaction, 1}(), Array{Metabolite, 1}(), Array{Gene, 1}())
end


"""
    getindex(model::CobraTools.Model, rxn::Reaction)

Get the index of rxn in model. Return -1 if not found.
"""
function Base.getindex(model::CobraTools.Model, rxn::Reaction)
    return model.reactions[rxn]
end

"""
    getindex(model::CobraTools.Model, met::Metabolite)

Get the index of metabolite in model. Return -1 if not found.
"""
function Base.getindex(model::CobraTools.Model, met::Metabolite)
    return model.metabolites[met]
end

"""
    getindex(model::CobraTools.Model, gene::Gene)

Get the index of gene in model. Return -1 if not found.
"""
function Base.getindex(model::CobraTools.Model, gene::Gene)
    return model.genes[gene]
end

"""
Pretty printing of model::CobraTools.Model.
"""
function Base.show(io::IO, ::MIME"text/plain", m::CobraTools.Model)
    println(io, "Constraint based model: ", m.id, "\n",
              "Number of reactions: ", length(m.reactions), "\n",
              "Number of metabolites: ", length(m.metabolites), "\n",
              "Number of genes: ", length(m.genes))
end


"""
    add!(model::CobraTools.Model, rxns::Union{Array{Reaction, 1}, Reaction})

Add rxn(s) to model if they are not already present (based on rxn.id)
"""
function add!(model::CobraTools.Model, rxns::Array{Reaction, 1})
    for rxn in rxns
        add!(model, rxn)
    end
end

function add!(model::CobraTools.Model, rxn::Reaction)
    if model[rxn] == -1
        push!(model.reactions, rxn)
    else
        @warn "$(rxn.id) already present in model."
    end
end

"""
    add!(model::CobraTools.Model, mets::Union{Array{Metabolite, 1}, Metabolite})

Add metabolite(s) to model if they are not already present (based on rxn.id)
"""
function add!(model::CobraTools.Model, mets::Array{Metabolite, 1})
    for met in mets
        add!(model, met)
    end
end

function add!(model::CobraTools.Model, met::Metabolite)
    if model[met] == -1
        push!(model.metabolites, met)
    else
        @warn "$(met.id) already present in model."
    end
end

"""
    add!(model::CobraTools.Model, genes::Union{Array{Gene, 1}, Gene})

Add gene(s) to model if they are not already present (based on gene.id)
"""
function add!(model::CobraTools.Model, genes::Array{Gene, 1})
    for gene in genes
        add!(model, gene)
    end
end

function add!(model::CobraTools.Model, gene::Gene)
    if model[gene] == -1
        push!(model.genes, gene)
    else
        @warn "$(gene.id) already present in model."
    end
end

"""
    rm!(model::CobraTools.Model, rxns::Union{Array{Reaction, 1}, Reaction})

Remove rxn(s) from model (based on ID).
"""
function rm!(model::CobraTools.Model, rxns::Union{Array{Reaction, 1}, Reaction})
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
end

"""
    rm!(model::CobraTools.Model, mets::Union{Array{Metabolite, 1}, Metabolite})

Remove met(s) from model (based on ID).
"""
function rm!(model::CobraTools.Model, mets::Union{Array{Metabolite, 1}, Metabolite})
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
end

"""
    rm!(model::CobraTools.Model, genes::Union{Array{Gene, 1}, Gene})

Remove gene(s) from model (based on ID).
"""
function rm!(model::CobraTools.Model, genes::Union{Array{Gene, 1}, Gene})
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
end

"""
    is_duplicate(model::CobraTools.Model, cmet::Metabolite)

Check if met already exists in model.metabolites but potentially has another id. 
First check if the id and formulas are the same. 
If not, check if the charges are the same.
If not, check if any of the annotations are the same.
"""
function is_duplicate(model::CobraTools.Model, cmet::Metabolite)
    _is_duplicate(model.metabolites, cmet)
end

"""
    is_duplicate(model::CobraTools.Model, crxn::Reaction)

Check if rxn already exists in model.reactions but potentially has another id. 
First checks if the ID already exists.
Then looks through the reaction equations and compares met.id's and stoichiometric coefficients.
If rxn has the same reaction equation as another reaction in rxns, the return true and the index of the match. 
"""
function is_duplicate(model::CobraTools.Model, crxn::Reaction)
    _is_duplicate(model.reactions, crxn)
end

"""
    is_duplicate(model::CobraTools.Model, cgene::Gene)

Check if gene already exists in model.genes but potentially has another id. 
First check if the ids are the same. 
If not, check if any of the annotations are the same.
"""
function is_duplicate(model::CobraTools.Model, cgene::Gene)
    _is_duplicate(model.genes, cgene)
end

"""
    fix_model!(model::CobraTools.Model)

Inspect metabolites and genes of model relative to the reactions.
Remove genes or metabolites that are not used in the reactions.
Add genes or metabolites that are not present in model.genes or model.metabolites but are used in model.reactions. 
"""
function fix_model!(model::CobraTools.Model)
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
        for gene_lists in rxn.grr
            for gene_list in gene_lists
                for gene in gene_list
                    if rxn_genes[gene] == -1
                        push!(rxn_genes, gene)
                    end
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
    !isempty(extra_mets) && rm!(model, [findfirst(model.reactions, x) for x in extra_mets])

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
