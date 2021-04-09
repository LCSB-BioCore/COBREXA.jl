module COBREXA

using Logging
using SparseArrays
using DelimitedFiles
using LinearAlgebra
using JuMP
using MAT
using Distributed
using DistributedData
using Downloads
using Requires
using JSON
using MacroTools
using Measurements
using Statistics
using Random
using Tulip # for LPs
using OSQP # for QPs, but it kinda sucks
using MacroTools # for DSL :)

import Base: findfirst, getindex, show
import Pkg
import SBML # conflict with Reaction struct name

include("banner.jl")
_print_banner()

# autoloading
const _inc(path...) = include(joinpath(path...))
const _inc_all(dir) = _inc.(joinpath.(dir, filter(fn -> endswith(fn, ".jl"), readdir(dir))))

_inc_all.(
    joinpath.(
        @__DIR__,
        [
            "types",
            "base",
            "io",
            joinpath("io", "show"),
            "sampling",
            "reconstruction",
            "analysis",
            "mods",
            "utils",
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

∅ = Metabolite("∅") # for exchange reactions
export ∅

end # module
