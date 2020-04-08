module ***REMOVED***

export LinearModel, addReaction, addReactions, solveLP, loadModel,
    nReactions, nMetabolites

include("modeling.jl")
include("solveLP.jl")
include("io/matReader.jl")
end
