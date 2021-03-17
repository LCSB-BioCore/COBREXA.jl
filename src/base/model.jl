"""
CobraModel struct of a constraint based metabolic model.

# Fields
````
id :: String
reactions :: Array{Reaction, 1}
metabolites :: Array{Metabolite, 1}
genes :: Array{Gene, 1}
````
"""
mutable struct CobraModel
    id::String
    reactions::Array{Reaction,1}
    metabolites::Array{Metabolite,1}
    genes::Array{Gene,1}
end

"""
Pretty printing of model::CobraModel.
"""
function Base.show(io::IO, ::MIME"text/plain", m::CobraModel)
    println(
        io,
        "Constraint based model: ",
        m.id,
        "\n",
        "Number of reactions: ",
        length(m.reactions),
        "\n",
        "Number of metabolites: ",
        length(m.metabolites),
        "\n",
        "Number of genes: ",
        length(m.genes),
    )
end

"""
CobraModel()

Empty model constructor.
"""
function CobraModel()
    CobraModel("blank", Array{Reaction,1}(), Array{Metabolite,1}(), Array{Gene,1}())
end

"""
    getindex(model::CobraModel, rxn::Reaction)

Get the index of `rxn` in `model`, based on reaction `id`. 
Return -1 if not found.

Typical usage: ind = model[rxn]
"""
function Base.getindex(model::CobraModel, rxn::Reaction)
    return model.reactions[rxn]
end

"""
    getindex(model::CobraModel, met::Metabolite)

Get the index of `met` in `model`, based on metabolite `id`. 
Return -1 if not found.

Typical usage: ind = model[met]
"""
function Base.getindex(model::CobraModel, met::Metabolite)
    return model.metabolites[met]
end

"""
    getindex(model::CobraModel, gene::Gene)

Get the index of `gene` in `model`, based on gene `id`. 
Return -1 if not found.

Typical usage: ind = model[gene]
"""
function Base.getindex(model::CobraModel, gene::Gene)
    return model.genes[gene]
end
