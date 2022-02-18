module COBREXA

using Distributed
using DistributedData
using JSON
using JuMP
using LinearAlgebra
using MAT
using MacroTools
using OrderedCollections
using Random
using Serialization
using SparseArrays
using StableRNGs
using Statistics

import Base: findfirst, getindex, show
import Pkg
import SBML # conflict with Reaction struct name

include("banner.jl")

# autoloading
const _inc(path...) = include(joinpath(path...))
const _inc_all(dir) = _inc.(joinpath.(dir, filter(fn -> endswith(fn, ".jl"), readdir(dir))))

_inc_all.(
    joinpath.(
        @__DIR__,
        [
            joinpath("base", "types", "abstract"),
            joinpath("base", "logging"),
            joinpath("base", "macros"),
            joinpath("base", "types"),
            "base",
            "io",
            joinpath("io", "show"),
            "reconstruction",
            joinpath("reconstruction", "modifications"),
            "analysis",
            joinpath("analysis", "modifications"),
            joinpath("analysis", "sampling"),
            joinpath("base", "utils"),
        ],
    ),
)

# export everything that isn't prefixed with _ (inspired by JuMP.jl, thanks!)
for sym in names(@__MODULE__, all = true)
    if sym in [Symbol(@__MODULE__), :eval, :include] || startswith(string(sym), ['_', '#'])
        continue
    end
    @eval export $sym
end

end # module
