
"""
Model fields define the basic information necessary to run analysis tools
"""
struct Model
    id :: String # model name
    S :: SparseMatrixCSC{Float64,Int64} # stoichiometric matrix
    b :: SparseVector{Float64,Int64} # mass balance rhs
    lbs :: Array{Float64, 1} # reaction lower bounds
    ubs :: Array{Float64, 1} # rxn upper bounds
    rxns :: Array{String, 1} # reactions
    mets :: Array{String, 1} # metabolites
    grrs :: Dict{String, Array{Array{String, 1}, 1}} # reaction -> [[genes]]
end

"""
Model()

Empty constructor.
"""
function Model()
    Model("blank",
            sparse(rand(0,0)), 
            sparse(rand(0)), 
            Float64[], 
            Float64[], 
            String[], 
            String[], 
            Dict{String, Array{Array{String, 1}, 1}}())
end

"""
Pretty printing of model.
"""
function Base.show(io::IO, m::Model)
    println(io, "Constraint based model: $(m.id)")
    println(io, "# rxns: $(length(m.rxns))")
    println(io, "# mets: $(length(m.mets))")
end