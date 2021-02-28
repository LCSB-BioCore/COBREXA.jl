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
using Tulip # for LPs
using OSQP # for QPs - only pFBA, not a very good LP solver

import Base: findfirst, getindex, show

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
include(joinpath("optimization_analysis", "basic_analysis.jl"))
# include("sampling_tools.jl")
# include("brenda_tools.jl")
# include("equilibrator_tools.jl")

∅ = Metabolite("∅") # for exchange reactions

export
    ∅,

    # base
    Model, Gene, Reaction, Metabolite, check_duplicate_annotations, get_atoms, check_same_formula, is_mass_balanced, check_duplicate_reaction,
    
    # construction
    ⟶, →, ←, ⟵, ↔, ⟷, add!, rm!, fix_model!,
    
    # io
    read_model, save_model,
    
    # optimization_analysis
    get_core_model, build_cbm, fba, map_fluxes, set_bound, pfba#, exchange_reactions, metabolite_fluxes

# Initialization functions
include("init_functions.jl")

end # module
