"""
StandardModel struct of a constraint based metabolic model.

# Fields
````
id :: String
reactions :: Vector{Reaction}
metabolites :: Vector{Metabolite}
genes :: Vector{Gene}
````
"""
mutable struct StandardModel <: MetabolicModel
    id::String
    reactions::Vector{Reaction}
    metabolites::Vector{Metabolite}
    genes::Vector{Gene}

    StandardModel(
        id = "",
        reactions::Vector{Reaction} = Reaction[],
        metabolites::Vector{Metabolite} = Metabolite[],
        genes::Vector{Gene} = Gene[],
    ) = new(id, reactions, metabolites, genes)
end

# MetabolicModel interface follows
reactions(model::StandardModel)::Vector{String} = [r.id for r in model.reactions]
n_reactions(model::StandardModel)::Int = length(model.reactions)

metabolites(model::StandardModel)::Vector{String} = [m.id for m in model.metabolites]
n_metabolites(model::StandardModel)::Int = length(model.metabolites)

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

balance(model::StandardModel)::SparseVec = spzeros(length(model.metabolites))

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
