"""
```
\\\\\\\\\\  // //     | COBREXA.jl  v$(COBREXA.COBREXA_VERSION)
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
using DocStringExtensions

import Base: findfirst, getindex, show
import Pkg
import SBML # conflict with Reaction struct name

const _PKG_ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
include_dependency(joinpath(_PKG_ROOT_DIR, "Project.toml"))

const COBREXA_VERSION =
    VersionNumber(Pkg.TOML.parsefile(joinpath(_PKG_ROOT_DIR, "Project.toml"))["version"])

end # module
