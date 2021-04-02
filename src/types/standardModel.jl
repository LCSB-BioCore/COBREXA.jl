"""
StandardModel struct of a constraint based metabolic model.

# Fields
````
id :: String
reactions :: Array{Reaction, 1}
metabolites :: Array{Metabolite, 1}
genes :: Array{Gene, 1}
````
"""
mutable struct StandardModel <: MetabolicModel
    id::String
    reactions::Array{Reaction,1}
    metabolites::Array{Metabolite,1}
    genes::Array{Gene,1}
end

"""
Pretty printing of model::StandardModel.
"""
function Base.show(io::IO, ::MIME"text/plain", m::StandardModel)
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
StandardModel()

Empty model constructor.
"""
function StandardModel()
    StandardModel("blank", Array{Reaction,1}(), Array{Metabolite,1}(), Array{Gene,1}())
end

"""
    getindex(model::StandardModel, rxn::Reaction)

Get the index of `rxn` in `model`, based on reaction `id`.
Return -1 if not found.

Typical usage: ind = model[rxn]
"""
function Base.getindex(model::StandardModel, rxn::Reaction)
    return model.reactions[rxn]
end

"""
    getindex(model::StandardModel, met::Metabolite)

Get the index of `met` in `model`, based on metabolite `id`.
Return -1 if not found.

Typical usage: ind = model[met]
"""
function Base.getindex(model::StandardModel, met::Metabolite)
    return model.metabolites[met]
end

"""
    getindex(model::StandardModel, gene::Gene)

Get the index of `gene` in `model`, based on gene `id`.
Return -1 if not found.

Typical usage: ind = model[gene]
"""
function Base.getindex(model::StandardModel, gene::Gene)
    return model.genes[gene]
end

# generic interface functions

function reactions(model::StandardModel)::Vector{String}
    [r.id for r in model.reactions]
end

function metabolites(model::StandardModel)::Vector{String}
    [m.id for m in model.metabolites]
end

function stoichiometry(model::StandardModel)::SparseMat
    S = SparseArrays.spzeros(length(model.metabolites), length(model.reactions))
    metids = metabolites(model)
    for (i, rxn) in enumerate(model.reactions) # column
        for (met, coeff) in rxn.metabolites
            j = findfirst(x -> x == met.id, metids) # row
            isnothing(j) ?
            (@error "S matrix construction error: $(met.id) not defined."; continue) :
            nothing
            S[j, i] = coeff
        end
    end
    return S
end

function bounds(model::StandardModel)::Tuple{SparseVec,SparseVec}
    ubs = [rxn.ub for rxn in model.reactions]
    lbs = [rxn.lb for rxn in model.reactions]
    return lbs, ubs
end

function balance(model::StandardModel)::SparseVec
    SparseArrays.spzeros(length(model.metabolites))
end

function objective(model::StandardModel)::SparseVec
    obj_arr = SparseArrays.spzeros(length(model.reactions))
    j = -1
    for (i, r) in enumerate(model.reactions)
        if r.objective_coefficient != 0.0
            j = i
            break
        end
    end
    if j != -1 # objective assigned, otherwise return array of 0s
        obj_arr[j] = 1.0
    end
    return obj_arr
end
