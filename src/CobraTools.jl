module CobraTools

using Requires
using MAT
using JSON
using SBML
using SparseArrays
using JuMP
using LinearAlgebra
using Measurements
using Statistics
using Random
using PyCall

# abstract types
abstract type ModelComponent end # for Reactions, Metabolites and Genes (all have IDs)

# definitions of structs
include("metabolite.jl")
include("gene.jl")
include("reaction.jl")
include("model.jl")

# Tools
include("io_tools.jl")
include("construction_overloading.jl")
include("brenda_tools.jl")
include("equilibrator_tools.jl")
# include("basic_analysis.jl")
# include("sampling_tools.jl")

# export
∅ = Metabolite("∅") # for exchange reactions
export ∅

export Model, add!, rm!, is_duplicate, fix_model! # from model
export Gene # from gene
export Metabolite, get_atoms # from metabolite
export Reaction, is_mass_balanced # from reaction
export ⟶, →, ←, ⟵, ↔, ⟷ # from construction_tools
export read_model, save_model # from io_tools
export build_cbm, fba, pfba, map_fluxes, set_bound, exchange_reactions, metabolite_fluxes


# Initialization functions
include("init_functions.jl")

end # module
