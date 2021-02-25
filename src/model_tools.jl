"""
Model struct of a constraint based metabolic model.
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
    Model("blank", Array{Reaction, 1}(), Array{Metabolite, 1}(), Array{Gene, 1}(), Dict{String, Array{Array{String, 1}, 1}}())
end


"""
S, b, upper_bounds, lower_bounds = get_core_model(model::Model)

Return stoichiometrix matrix (S), mass balance right hand side (b), upper and lower bounds of constraint based model.
That is, S*v=b with lower_bounds ≤ v ≤ upper_bounds where v is the flux vector. This is useful if you want to construct
your own optimization problem.
"""
function get_core_model(model::Model)
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
Pretty printing of model::Model.
"""
function Base.show(io::IO, m::Model)
    println(io, "Constraint based model: ", m.id)
    println(io, "Number of reactions: ", length(m.rxns))
    println(io, "Number of metabolites: ", length(m.mets))
    println(io, "Number of genes: ", length(m.genes))
end
