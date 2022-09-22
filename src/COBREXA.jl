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

import Pkg
import FromFile: @from

const _PKG_ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
include_dependency(joinpath(_PKG_ROOT_DIR, "Project.toml"))

const COBREXA_VERSION =
    VersionNumber(Pkg.TOML.parsefile(joinpath(_PKG_ROOT_DIR, "Project.toml"))["version"])

#=
Organize COBREXA into modules to make discovery of functions and data structures
easier. Load order matters inside each module. 
=#

"""
Common types, utilities and constants used by COBREXA.
"""
module Common # needed a synonym for Base :) 
    
    include("dependencies.jl") # TODO: fix this creates duplicate names, use FromFile

    include(joinpath("base", "types", "abstract", "Maybe.jl"))
    include(joinpath("base", "types", "abstract", "MetabolicModel.jl"))

    include(joinpath("base", "logging", "log.jl"))
    include(joinpath("base", "macros", "change_bounds.jl"))
    include(joinpath("base", "macros", "is_xxx_reaction.jl"))
    include(joinpath("base", "macros", "modell_wrapper.jl"))
    include(joinpath("base", "macros", "remove_item.jl"))
    include(joinpath("base", "macros", "serialized.jl"))
 
    include(joinpath("base", "types", "MetabolicModel.jl"))
    include(joinpath("base", "types", "CoreModel.jl"))
    include(joinpath("base", "types", "CoreModelCoupled.jl"))
    include(joinpath("base", "types", "Gene.jl"))
    include(joinpath("base", "types", "Reaction.jl"))
    include(joinpath("base", "types", "Metabolite.jl"))
    include(joinpath("base", "types", "StandardModel.jl"))

    include(joinpath("base", "types", "FluxSummary.jl"))
    include(joinpath("base", "types", "FluxVariabilitySummary.jl"))
    
    include(joinpath("base", "types", "HDF5Model.jl"))
    include(joinpath("base", "types", "JSONModel.jl"))
    include(joinpath("base", "types", "MATModel.jl"))
    include(joinpath("base", "types", "SBMLModel.jl"))


end
using .Common

"""
Input/output functions, as well as pretty printing.
"""
module IO
    include("dependencies.jl") # TODO: fix this creates duplicate names, use FromFile
    using ..Common

    # IO functions
    include(joinpath("io", "json.jl"))

end
using .IO

# export everything that isn't prefixed with _ (inspired by JuMP.jl, thanks!)
for sym in names(@__MODULE__, all = true)
    if sym in [Symbol(@__MODULE__), :eval, :include] || startswith(string(sym), ['_', '#'])
        continue
    end
    @eval export $sym
end

end # module
