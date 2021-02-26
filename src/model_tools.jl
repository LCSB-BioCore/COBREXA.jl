"""
Model struct of a constraint based metabolic model.

# Fields
````
id :: String
rxns :: Array{Reaction, 1}
mets :: Array{Metabolite, 1}
genes :: Array{Gene, 1}
grrs :: Dict{String, Array{Array{String, 1}, 1}}
````
"""
struct Model
    id :: String # model name
    rxns :: Array{Reaction, 1} # reaction metadata
    mets :: Array{Metabolite, 1}  # metabolite metadata
    genes :: Array{Gene, 1} # gene metadata
    grrs :: Dict{String, Array{Array{String, 1}, 1}} # reaction -> [[gene & gene...] or [gene & gene...] or [gene & gene...]] 
end

"""
Model()

Empty model constructor.
"""
function Model()
    CobraTools.Model("blank", Array{Reaction, 1}(), Array{Metabolite, 1}(), Array{Gene, 1}(), Dict{String, Array{Array{String, 1}, 1}}())
end

"""
index = getindex(model::CobraTools.Model, rxn::Reaction)

Get the index of rxn in model. Return -1 if not found.
"""
function Base.getindex(model::CobraTools.Model, rxn::Reaction)
    return model.rxns[rxn]
end

"""
index = getindex(model::CobraTools.Model, met::Metabolite)

Get the index of metabolite in model. Return -1 if not found.
"""
function Base.getindex(model::CobraTools.Model, met::Metabolite)
    return model.mets[met]
end

"""
index = getindex(model::CobraTools.Model, gene::Gene)

Get the index of gene in model. Return -1 if not found.
"""
function Base.getindex(model::CobraTools.Model, gene::Gene)
    return model.genes[gene]
end

"""
S, b, upper_bounds, lower_bounds = get_core_model(model::CobraTools.Model)

Return stoichiometrix matrix (S), mass balance right hand side (b), upper and lower bounds of constraint based model.
That is, S*v=b with lower_bounds ≤ v ≤ upper_bounds where v is the flux vector. This is useful if you want to construct
your own optimization problem.
"""
function get_core_model(model::CobraTools.Model)
    ubs = [rxn.ub for rxn in model.rxns]
    lbs = [rxn.lb for rxn in model.rxns]
    
    b = spzeros(length(model.mets))
    S = spzeros(length(model.mets), length(model.rxns))

    metids = [met.id for met in model.mets] # need indices for S matrix construction
    for (i, rxn) in enumerate(model.rxns) # column
        for (met, coeff) in rxn.metabolites
            j = findfirst(x -> x == met.id, metids) # row
            isnothing(j) ? (@error "S matrix construction error: $(met.id) not defined."; continue) : nothing
            S[j, i] = coeff
        end
    end
    return S, b, ubs, lbs
end

"""
Pretty printing of model::CobraTools.Model.
"""
function Base.show(io::IO, ::MIME"text/plain", m::CobraTools.Model)
    println(io, "Constraint based model: ", m.id, "\n",
              "Number of reactions: ", length(m.rxns), "\n",
              "Number of metabolites: ", length(m.mets), "\n",
              "Number of genes: ", length(m.genes))
end

