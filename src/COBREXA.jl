module COBREXA

using Logging
using SparseArrays
using DelimitedFiles
using LinearAlgebra
using JuMP
using MAT
using Distributed
using DistributedData
using Requires
using JSON
using Measurements
using Statistics
using Random
using PyCall
using Tulip # for LPs
using OSQP # for QPs, but it kinda sucks

import Base: findfirst, getindex, show
import Pkg
import SBML # conflict with Reaction struct name

include("banner.jl")
_printBanner()

# autoloading
const _inc(path...) = include(joinpath(path...))
const _inc_all(dir) = _inc.(joinpath.(dir, filter(fn -> endswith(fn, ".jl"), readdir(dir))))
_inc_all.(joinpath.(@__DIR__, ["types", "base","analysis"]))#"io", "reconstruction", "analysis"]))

# export everything that isn't prefixed with _ (inspired by JuMP.jl, thanks!)
for sym in names(@__MODULE__, all = true)
    if sym in [Symbol(@__MODULE__), :eval, :include] || startswith(string(sym), ['_', '#'])
        continue
    end
    @eval export $sym
end

end # module
