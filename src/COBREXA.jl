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
_printBanner()

const _inc(path...) = include(joinpath(path...))

_inc.("types", ["linearModel.jl", "reactionStatus.jl"])

_inc.("base", ["solver.jl", "utilities.jl"])
_inc.("io", ["reader.jl", "writer.jl", "sbml.jl"])
_inc.("reconstruction", ["coupling.jl", "modeling.jl"])
_inc.("analysis", ["fba.jl", "fva.jl"])


# export everything that isn't prefixed with _ (inspired by JuMP.jl, thanks!)
for sym in names(@__MODULE__, all = true)
    if sym in [Symbol(@__MODULE__), :eval, :include] || startswith(string(sym), "_")
        continue
    end
    @eval export $sym
end

end
