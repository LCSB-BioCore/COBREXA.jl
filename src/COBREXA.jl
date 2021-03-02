module COBREXA

using Logging
using SparseArrays
using DelimitedFiles
using LinearAlgebra
using JuMP
using MAT
using Distributed
using SBML
using DistributedData
import Pkg

include("banner.jl")
printBanner()

include(joinpath("header", "loadSource.jl"))
include(joinpath("header", "types.jl"))

loadSource(["base", "io", "reconstruction", "analysis"], @__DIR__)

# export functions
export speye, LinearModel, nReactions, nMetabolites, nCouplingConstraints,
       addReaction, addReactions, removeReactions, changeBounds!,
       addCouplingConstraints!, addCouplingConstraints,
       removeCouplingConstraints!, removeCouplingConstraints,
       changeCouplingBounds!, changeCouplingBounds,
       verifyConsistency, findExchangeReactions, findExchangeMetabolites,
       solveLP, loadModel, loadSBMLModel, writeModel,
       fluxBalanceAnalysis, fluxVariabilityAnalysis, parFVA,
       convertToExportable

end
