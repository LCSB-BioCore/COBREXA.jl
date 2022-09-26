"""
```
\\\\\\\\\\  // //     | COBREXA.jl  v$(COBREXA.version)
 \\\\ \\\\// //      |
  \\\\ \\/ //       | COnstraint-Based Reconstruction
   \\\\  //        | and EXascale Analysis in Julia
   //  \\\\        |
  // /\\ \\\\       | See documentation and examples at:
 // //\\\\ \\\\      | https://lcsb-biocore.github.io/COBREXA.jl
// //  \\\\\\\\\\     |
```

To start up quickly, install your favorite optimizer, load a metabolic model in
a format such as SBML or JSON, and run a metabolic analysis such as the flux
balance analysis:
```
import Pkg; Pkg.add("GLPK")
using COBREXA, GLPK
model = load_model("e_coli_core.xml")
x = flux_balance_analysis_dict(model, GLPK.Optimizer)
flux_summary(x)
```

A complete overview of the functionality can be found in the documentation.
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

# versioning tools
const _PKG_ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
include_dependency(joinpath(_PKG_ROOT_DIR, "Project.toml"))

const version =
    VersionNumber(Pkg.TOML.parsefile(joinpath(_PKG_ROOT_DIR, "Project.toml"))["version"])

module ModuleTools
macro inc(path...)
    esc(:(include(joinpath(@__DIR__, $(joinpath(String.(path)...) * ".jl")))))
end

macro inc_dir(path...)
    dir = joinpath(@__DIR__, String.(path)...)
    files = filter(endswith(".jl"), readdir(dir; join = true))
    esc(Expr(:block, (:(include($f)) for f in files)...))
end

macro dse()
    :(using DocStringExtensions)
end

# export everything from the local namespace that seems exportable
# (inspired by JuMP.jl, thanks!)
macro export_locals()
    quote
        for sym in names(@__MODULE__; all = true, imported = true)
            sym in [Symbol(@__MODULE__), :eval, :include] && continue
            startswith(string(sym), ['_', '#']) && continue
            sym == :Internal && continue
            @eval export $(Expr(:$, :sym))
        end
    end
end

@export_locals
end

# start loading the individual modules
module Internal
using ..ModuleTools
@inc macros
end

using .ModuleTools

@inc log
@inc types
@inc io

#= TODOs here

_inc_all.(
    joinpath.(
        @__DIR__,
        [
            joinpath("types", "abstract"),
            joinpath("logging"),
            joinpath("macros"),
            joinpath("types"),
            joinpath("types", "wrappers"),
            joinpath("ontology"),
            "base",
            "io",
            joinpath("io", "show"),
            "reconstruction",
            joinpath("reconstruction", "modifications"),
            "analysis",
            joinpath("analysis", "modifications"),
            joinpath("analysis", "sampling"),
            joinpath("utils"),
        ],
    ),
)
=#

end # module
