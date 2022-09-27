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

import Pkg

# versioning tools
const _PKG_ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
include_dependency(joinpath(_PKG_ROOT_DIR, "Project.toml"))

const version =
    VersionNumber(Pkg.TOML.parsefile(joinpath(_PKG_ROOT_DIR, "Project.toml"))["version"])

# bootstrap the module machinery
include("modules.jl")
using .ModuleTools

# load various internal helpers first
@inc internal
@inc log

# start loading individual user-facing modules
@inc types

@inc io
@inc solver
@inc reconstruction

# TODO: this needs to be assimilated to actual modules
module Rest
using ..ModuleTools
@dse

using ..Internal
using ..Types
using ..Accessors
using ..IO

#"analysis",
#joinpath("analysis", "modifications"),
#joinpath("analysis", "sampling"),
#joinpath("utils"),
end

end # module
