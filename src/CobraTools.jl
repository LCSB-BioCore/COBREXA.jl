module CobraTools

# model reading and writing
using JSON
using MATLAB 
using PyCall # NB: need to install libsbml

# Model analysis
using SparseArrays
using JuMP

# set the output flag (turns on/true info and warning messages)
verbose = true

include("cobra.jl")
include("parsemodels.jl")
include("analysis.jl")
include("rxn_construction.jl")

# cobra
export Reaction, Metabolite, Gene

# rxn_construction
∅ = Metabolite("∅") # for exchange reactions
export ∅, ⟶, →, ←, ⟵, ↔, ⟷

end # module
