module COBREXA

using Logging
using SparseArrays
using DelimitedFiles
using LinearAlgebra
using JuMP
using MAT
using Distributed
using SBML
using DistributedData
import Pkg

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


loadSource(["base", "io", "reconstruction", "analysis"], @__DIR__)

# export functions
export speye, LinearModel, nReactions, nMetabolites, nCouplingConstraints,
       addReaction, addReactions, removeReactions, changeBounds!,
       addCouplingConstraints!, addCouplingConstraints,
       removeCouplingConstraints!, removeCouplingConstraints,
       changeCouplingBounds!, changeCouplingBounds,
       verifyConsistency, findExchangeReactions, findExchangeMetabolites,
       solveLP, loadModel, loadSBMLModel, writeModel,
       fluxBalanceAnalysis, fluxVariabilityAnalysis, parFVA,
       convertToExportable

end


##########################################################
using Requires
using JSON
using Measurements
using Statistics, Random
using PyCall
using JuMP, Tulip, OSQP # OSQP sucks for LPs

import Base: findfirst, getindex, show

# definitions of structs
include(joinpath("base", "metabolite.jl"))
include(joinpath("base", "gene.jl"))
include(joinpath("base", "reaction.jl"))
include(joinpath("base", "model.jl"))

# Tools
include(joinpath("io", "io_tools.jl"))
include(joinpath("construction", "construction_overloading.jl"))
include(joinpath("construction", "model_manipulations.jl"))
include(joinpath("optimization_analysis", "basic_analysis.jl"))
include(joinpath("sampling", "sampling_tools.jl"))
include(joinpath("external", "brenda_tools.jl"))
include(joinpath("external", "equilibrator_tools.jl"))

∅ = Metabolite("∅") # for exchange reactions

export ∅,

    # base
    Model,
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
    achr,

    # external
    parse_brenda,
    map_gibbs_rxns,
    map_gibbs_external,
    map_gibbs_internal

# Initialization functions
include("init_functions.jl")