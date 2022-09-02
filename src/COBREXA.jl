"""
```
\\\\\\\\\\  // //     | COBREXA.jl  v$(VersionNumber(Pkg.TOML.parsefile(joinpath(normpath(joinpath(@__DIR__, "..")), "Project.toml"))["version"]))
 \\\\ \\\\// //      |
  \\\\ \\/ //       | COnstraint-Based Reconstruction
   \\\\  //        | and EXascale Analysis in Julia
   //  \\\\        |
  // /\\ \\\\       | See documentation and examples at:
 // //\\\\ \\\\      | https://lcsb-biocore.github.io/COBREXA.jl
// //  \\\\\\\\\\     |
```
To perform flux balance analysis using COBREXA, try:
```
using COBREXA, Clarabel

model = load_model(StandardModel, "path_to_model")

fba_dict = flux_balance_analysis_dict(model, Clarabel.Optimizer)

flux_summary(flux_dict)
```
More examples can be found in the documentation.
"""
module COBREXA

using Distributed
using DistributedData
using HDF5
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

const _PKG_ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
include_dependency(joinpath(_PKG_ROOT_DIR, "Project.toml"))

const COBREXA_VERSION =
    VersionNumber(Pkg.TOML.parsefile(joinpath(_PKG_ROOT_DIR, "Project.toml"))["version"])

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
            joinpath("base", "types", "wrappers"),
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
