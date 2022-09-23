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
import SBML # conflict with Reaction struct name 
import Pkg

include_dependency(joinpath(normpath(joinpath(@__DIR__, "..")), "Project.toml"))
const COBREXA_VERSION =
    VersionNumber(Pkg.TOML.parsefile(joinpath(normpath(joinpath(@__DIR__, "..")), "Project.toml"))["version"])

#=
Organize COBREXA into modules to make discovery of functions and data structures
easier. Load order matters inside each module. 
=#

module Common
using ..COBREXA.DocStringExtensions

include(joinpath("base", "types", "abstract", "Maybe.jl"))

end
using .Common

module Macros
using ..COBREXA.DocStringExtensions

include(joinpath("base", "logging", "log.jl"))
include(joinpath("base", "macros", "model_wrapper.jl"))
include(joinpath("base", "macros", "is_xxx_reaction.jl"))
include(joinpath("base", "macros", "change_bounds.jl"))
include(joinpath("base", "macros", "remove_item.jl"))
end
using .Macros

# macros
# include(joinpath("base", "macros", "serialized.jl"))


# constants and definitions
include(joinpath("base", "constants.jl"))
include(joinpath("base", "SBOTerms.jl")) # must come before identifiers.jl
include(joinpath("base", "identifiers.jl"))

"""
Common types, utilities and constants used by COBREXA.
"""
module ModelTypes
using ..Common: Maybe
using ..Macros: @_io_log, @_inherit_model_methods, @_inherit_model_methods_fn
using ..COBREXA: _constants
using ..COBREXA.JSON,
    ..COBREXA.SparseArrays,
    ..COBREXA.DocStringExtensions,
    ..COBREXA.MAT,
    ..COBREXA.HDF5,
    ..COBREXA.SBML,
    ..COBREXA.OrderedCollections

# abstract models
include(joinpath("base", "types", "abstract", "MetabolicModel.jl"))
include(joinpath("base", "types", "MetabolicModel.jl"))

# concrete external models
include(joinpath("base", "types", "JSONModel.jl"))
include(joinpath("base", "types", "HDF5Model.jl"))
include(joinpath("base", "types", "MATModel.jl"))
include(joinpath("base", "types", "SBMLModel.jl"))

# concrete internal models and types
include(joinpath("base", "types", "CoreModel.jl"))
include(joinpath("base", "types", "CoreModelCoupled.jl"))

include(joinpath("base", "types", "Gene.jl"))
include(joinpath("base", "types", "Reaction.jl"))
include(joinpath("base", "types", "Metabolite.jl"))
include(joinpath("base", "types", "StandardModel.jl"))

include(joinpath("base", "types", "Serialized.jl"))

include(joinpath("base", "types", "ReactionStatus.jl"))

include(joinpath("base", "types", "ModelWrapper.jl"))

# enzyme constrained models 
include(joinpath("base", "types", "Isozyme.jl"))
include(joinpath("base", "types", "wrappers", "GeckoModel.jl"))
include(joinpath("base", "types", "wrappers", "SMomentModel.jl"))

# io types
include(joinpath("base", "types", "FluxSummary.jl"))
include(joinpath("base", "types", "FluxVariabilitySummary.jl"))
export CoreModel
end
using .ModelTypes
export CoreModel

"""
Input/output functions, as well as pretty printing.
"""
module InputOutput # can't use IO
using ..ModelTypes:
    JSONModel,
    SBMLModel,
    MATModel,
    MetabolicModel,
    HDF5Model,
    FluxSummary,
    FluxVariabilitySummary,
    Gene,
    Metabolite,
    Reaction,
    Serialized
using ..Macros: @_io_log
using ..COBREXA.JSON,
    ..COBREXA.MAT, ..COBREXA.SBML, ..COBREXA.HDF5, ..COBREXA.DocStringExtensions

# IO functions
include(joinpath("io", "json.jl"))
include(joinpath("io", "mat.jl"))
include(joinpath("io", "sbml.jl"))
include(joinpath("io", "h5.jl"))
include(joinpath("io", "io.jl"))

# pretty printing 
include(joinpath("io", "show", "pretty_printing.jl"))
include(joinpath("io", "show", "FluxSummary.jl"))
include(joinpath("io", "show", "FluxVariabilitySummary.jl"))
include(joinpath("io", "show", "Gene.jl"))
include(joinpath("io", "show", "MetabolicModel.jl"))
include(joinpath("io", "show", "Metabolite.jl"))
include(joinpath("io", "show", "Reaction.jl"))
include(joinpath("io", "show", "Serialized.jl"))

end
using .InputOutput

module Utils
using ..ModelTypes:
    Gene,
    Metabolite,
    Reaction,
    MetaboliteFormula,
    CoreModel,
    CoreModelCoupled,
    GeckoModel,
    SMomentModel,
    MetabolicModel,
    GeneAssociation,
    Isozyme,
    MetabolicModel,
    StandardModel,
    Serialized

using ..Macros: @_models_log, @_is_reaction_fn

using ..COBREXA.SBML,
    ..COBREXA.HDF5,
    ..COBREXA.OrderedCollections,
    ..COBREXA.Serialization,
    ..COBREXA.DocStringExtensions,
    ..COBREXA.SparseArrays

include(joinpath("base", "utils", "Annotation.jl"))
include(joinpath("base", "utils", "bounds.jl"))
include(joinpath("base", "utils", "chemical_formulas.jl"))
include(joinpath("base", "utils", "CoreModel.jl"))
include(joinpath("base", "utils", "enzymes.jl"))
include(joinpath("base", "utils", "fluxes.jl"))
include(joinpath("base", "utils", "gecko.jl"))
include(joinpath("base", "utils", "gene_associations.jl"))
include(joinpath("base", "utils", "guesskey.jl"))
include(joinpath("base", "utils", "HDF5Model.jl"))
include(joinpath("base", "utils", "looks_like.jl"))
include(joinpath("base", "utils", "Reaction.jl"))
include(joinpath("base", "utils", "Serialized.jl"))
include(joinpath("base", "utils", "smoment.jl"))
include(joinpath("base", "utils", "StandardModel.jl"))
end
using .Utils

module Analysis
using ..ModelTypes:
    Gene,
    Metabolite,
    Reaction,
    CoreModel,
    CoreModelCoupled,
    GeckoModel,
    SMomentModel,
    MetabolicModel,
    Isozyme,
    MetabolicModel,
    StandardModel,
    Serialized,
    SparseMat

using ..Macros: @_models_log

using ..Common: Maybe

using ..COBREXA.DocStringExtensions,
    ..COBREXA.JuMP,
    ..COBREXA.Distributed,
    ..COBREXA.DistributedData

include(joinpath("base", "solver.jl"))
include(joinpath("analysis", "screening.jl"))

include(joinpath("analysis", "flux_balance_analysis.jl"))
include(joinpath("analysis", "flux_variability_analysis.jl"))
include(joinpath("analysis", "minimize_metabolic_adjustment.jl"))
include(joinpath("analysis", "parsimonious_flux_balance_analysis.jl"))
include(joinpath("analysis", "max_min_driving_force.jl"))
include(joinpath("analysis", "gecko.jl"))
include(joinpath("analysis", "smoment.jl"))

include(joinpath("analysis", "sampling", "warmup_variability.jl"))
include(joinpath("analysis", "sampling", "affine_hit_and_run.jl"))

module Modifications # optimization modifications
using ....ModelTypes:
    MetabolicModel

using ....COBREXA.DocStringExtensions,
    ....COBREXA.JuMP
    
include(joinpath("analysis", "modifications", "generic.jl"))
include(joinpath("analysis", "modifications", "crowding.jl"))
include(joinpath("analysis", "modifications", "knockout.jl"))
include(joinpath("analysis", "modifications", "loopless.jl"))
include(joinpath("analysis", "modifications", "moment.jl")) # TODO remove and deprecate
include(joinpath("analysis", "modifications", "optimizer.jl"))
end
using .Modifications
end
using .Analysis
using .Analysis.Modifications

module Reconstruction

using ..ModelTypes:
    Gene,
    Metabolite,
    Reaction,
    CoreModel,
    CoreModelCoupled,
    MetabolicModel,
    MetabolicModel,
    StandardModel,
    SparseMat,
    StringVecType,
    VecType,
    MatType

using ..Macros: @_change_bounds_fn, @_remove_fn

using ..Common: Maybe

using ..COBREXA.DocStringExtensions

include(joinpath("reconstruction", "CoreModel.jl")) # borked

end
using .Reconstruction
using .Reconstruction.Modifications


# export everything that isn't prefixed with _ (inspired by JuMP.jl, thanks!)
for sym in names(@__MODULE__, all = true)
    if sym in [Symbol(@__MODULE__), :eval, :include] || startswith(string(sym), ['_', '#'])
        continue
    end
    @eval export $sym
end

end # module
