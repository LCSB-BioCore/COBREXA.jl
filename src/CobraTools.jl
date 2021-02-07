module CobraTools

# model reading and writing
using JSON
using MATLAB 
using PyCall # NB: need to install libsbml

# Model analysis
using SparseArrays
using JuMP

# Global options for package
mutable struct CobraToolsOptions
    verbose :: Bool
end
cto = CobraToolsOptions(true)

"""
Reduce verbosity (@info and @warn are suppressed)

This is mostly useful for unit testing to suppress all the model reading and writing that occurs therein.
"""
function setverbose(verbose)
    if verbose
        cto.verbose = true
    else
        cto.verbose = false
    end
end

include("cobra.jl")
export Reaction, Metabolite, Gene

include("parsemodels.jl")

include("analysis.jl")

include("rxn_construction.jl")
∅ = Metabolite("∅") # for exchange reactions
export ∅, ⟶, →, ←, ⟵, ↔, ⟷


end # module
