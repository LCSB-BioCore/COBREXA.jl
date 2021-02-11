module CobraTools

# IO of models and data
using JSON
using MATLAB 
using JLD

# Model analysis
using SparseArrays
using JuMP
using LinearAlgebra
# Find a way to only import packages the user actually has...?
using Gurobi
using Tulip
using GLPK
using Ipopt

# Gibbs
using Measurements
using Statistics
using PyCall # for Equilibrator - ensure that it is installed

include("global_cobratools.jl")

include("cobra.jl")
export Reaction, Metabolite, Gene

include("parse_models.jl")

include("rxn_tools.jl")
∅ = Metabolite("∅") # for exchange reactions
export ∅, ⟶, →, ←, ⟵, ↔, ⟷

include("met_tools.jl")

include("basic_analysis.jl")
# export Solution

include("gibbs_tools.jl")
include("name_space.jl")


# Init function
function __init__()
    py"""
    from equilibrator_api import ComponentContribution, Q_
    cc = ComponentContribution()

    def pygetdg0(f, ph, ionic):
        rxn = cc.parse_reaction_formula(f)
        isbal = rxn.is_balanced()
        cc.p_h = Q_(ph)  # set pH
        cc.ionic_strength = Q_(ionic)  # set I
        return isbal, cc.standard_dg(rxn).value.magnitude, cc.standard_dg(rxn).error.magnitude
    
    def pygetdgprime(f, ph, ionic):
        rxn = cc.parse_reaction_formula(f)
        isbal = rxn.is_balanced()
        cc.p_h = Q_(ph)  # set pH
        cc.ionic_strength = Q_(ionic)  # set I
        return isbal, cc.standard_dg_prime(rxn).value.magnitude, cc.standard_dg_prime(rxn).error.magnitude

    def pygetdgprimephys(f, ph, ionic):
        rxn = cc.parse_reaction_formula(f)
        isbal = rxn.is_balanced()
        cc.p_h = Q_(ph)  # set pH
        cc.ionic_strength = Q_(ionic)  # set I
        return isbal, cc.physiological_dg_prime(rxn).value.magnitude, cc.physiological_dg_prime(rxn).error.magnitude
    """


end

end # module
