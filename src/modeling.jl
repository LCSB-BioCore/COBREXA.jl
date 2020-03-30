mutable struct CobraLP
    S       ::AbstractMatrix
    b       ::Array{Float64,1}
    c       ::Array{Float64,1}
    lb      ::Array{Float64,1}
    ub      ::Array{Float64,1}
    rxns    ::Array{String,1}
    mets    ::Array{String,1}

    function CobraLP(
    S       ::AbstractMatrix,
    b       ::Array{Float64,1},
    c       ::Array{Float64,1},
    lb      ::Array{Float64,1},
    ub      ::Array{Float64,1},
    rxns    ::Array{String,1},
    mets    ::Array{String,1})

    new(S,b,c,lb,ub,rxns,mets)
    end
end
