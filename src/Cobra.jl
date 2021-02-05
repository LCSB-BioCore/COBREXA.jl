
"""
Model fields define the basic information necessary to run analysis tools
"""
struct Model
    S :: SparseMatrixCSC{Float64,Int64} # stoichiometric matrix
    b :: SparseVector{Float64,Int64} # mass balance rhs
    lb :: Array{Float64, 1} # reaction lower bounds
    ub :: Array{Float64, 1} # rxn upper bounds
    rxns :: Array{String, 1} # reactions
    mets :: Array{String, 1} # metabolites
    grrs # gene-reaction-rules
end

