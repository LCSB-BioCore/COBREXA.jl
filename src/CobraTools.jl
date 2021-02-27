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
abstract type ModelComponent end # for Reactions, Metabolites and Genes. All ModelComponents have an `id` field.

# definitions of structs
include(joinpath("base", "metabolite.jl"))
include(joinpath("base", "gene.jl"))
include(joinpath("base", "reaction.jl"))
include(joinpath("base", "model.jl"))

# Tools
include(joinpath("io", "io_tools.jl"))
include(joinpath("construction", "construction_overloading.jl"))
include(joinpath("construction", "model_manipulations.jl"))
# include("brenda_tools.jl")
# include("equilibrator_tools.jl")
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
export build_cbm, fba, pfba, map_fluxes, set_bound, exchange_reactions, metabolite_fluxes, get_core_model

# Initialization functions
include("init_functions.jl")

end # module
