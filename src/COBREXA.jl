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

# import src files
include(joinpath("header", "header.jl"))
include(joinpath("header", "types.jl"))

const PKG_ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
include_dependency(joinpath(PKG_ROOT_DIR, "Project.toml"))
const COBREXA_VERSION = VersionNumber(Pkg.TOML.parsefile(joinpath(PKG_ROOT_DIR, "Project.toml"))["version"])

c = Base.text_colors
tx = c[:normal] # text
b = c[:bold] * c[:blue]
r = c[:bold] * c[:red]
g = c[:bold] * c[:green]
m = c[:bold] * c[:magenta]
banner = "
   ____ ___  ____  ____  _____$(g)__$(tx)  $(r)__$(tx)    _     |
  / ___/ _ \\| __ )|  _ \\| ____$(g)\\ \\$(tx)$(r)/ /$(tx)   / \\    | Constraint-Based Reconstruction
 | |  | | | |  _ \\| |_) |  _|  $(g)\\$(tx)  $(r)/$(tx)   / _ \\   | and EXascale Analysis in Julia
 | |__| |_| | |_) |  _ <| |___ $(m)/$(tx)  $(b)\\$(tx)  / ___ \\  |
  \\____\\___/|____/|_| \\_\\_____$(m)/_/$(tx)$(b)\\_\\$(tx)/_/   \\_\\ | Version: v$(COBREXA_VERSION)
                                              |
"

print(banner)

# order issue
include("base/gene.jl")
include("base/metabolite.jl")
include("base/reaction.jl")
include("base/model.jl")
include("base/solver.jl")
include("base/utilities.jl")

loadSource(["io", "reconstruction", "analysis", "sampling"], @__DIR__)

# export functions
∅ = Metabolite("∅") # for exchange reactions
export speye, LinearModel, nReactions, nMetabolites, nCouplingConstraints,
       addReaction, addReactions, removeReactions, changeBounds!,
       addCouplingConstraints!, addCouplingConstraints,
       removeCouplingConstraints!, removeCouplingConstraints,
       changeCouplingBounds!, changeCouplingBounds,
       verifyConsistency, findExchangeReactions, findExchangeMetabolites,
       solveLP, loadModel, loadSBMLModel, writeModel,
       fluxBalanceAnalysis, fluxVariabilityAnalysis, parFVA,
       convertToExportable,
    
       # from CobraTools 
    ∅,
    # base
    CobraModel,
    Gene,
    Reaction,
    Metabolite,
    check_duplicate_annotations,
    get_atoms,
    check_same_formula,
    is_mass_balanced,
    check_duplicate_reaction,

    # construction
    ⟶,
    →,
    ←,
    ⟵,
    ↔,
    ⟷,
    add!,
    rm!,
    fix_model!,

    # io
    read_model,
    save_model,

    # optimization_analysis
    get_core_model,
    build_cbm,
    fba,
    map_fluxes,
    set_bound,
    pfba,
    atom_exchange,
    exchange_reactions,
    metabolite_fluxes,
    fva,

    # sampling
    hit_and_run,
    test_samples,
    achr

end # module