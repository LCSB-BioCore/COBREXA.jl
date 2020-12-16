module COBREXA

using Logging
using SparseArrays
using DelimitedFiles
using LinearAlgebra
using JuMP
using MAT
using Distributed
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
       solveLP, loadModel, fluxBalanceAnalysis, fluxVariabilityAnalysis, parFVA,
       writeModel, convertToExportable, createParPool

end
