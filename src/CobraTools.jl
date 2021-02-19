module CobraTools

# IO of models and data
using JSON
using MATLAB
using SBML

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

# Sampling
using Random

include("cobra_base.jl")
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

include("sampling.jl")


# Init function - build Gibbs calling functions
function __init__()
    py"""
    from equilibrator_api import ComponentContribution, Q_
    cc = ComponentContribution()

    def pygetdg0(fs, ph, ionic):
        cc.p_h = Q_(ph)  # set pH
        cc.ionic_strength = Q_(ionic)  # set I
        bals = []
        mags = []
        errs = []
        for f in fs:
            try:
                rxn = cc.parse_reaction_formula(f)
                isbal = rxn.is_balanced()
                bals += [isbal]
                mags += [cc.standard_dg(rxn).value.magnitude]
                errs += [cc.standard_dg(rxn).error.magnitude]
            except:
                bals += [False]
                mags += [0.0]
                errs += [0.0]
                
        return bals, mags, errs
    
    def pygetdgprime(fs, ph, ionic):
        cc.p_h = Q_(ph)  # set pH
        cc.ionic_strength = Q_(ionic)  # set I
        bals = []
        mags = []
        errs = []
        for f in fs:
            try:
                rxn = cc.parse_reaction_formula(f)
                isbal = rxn.is_balanced()
                bals += [isbal]
                mags += [cc.standard_dg_prime(rxn).value.magnitude]
                errs += [cc.standard_dg_prime(rxn).error.magnitude]
            except:
                bals += [False]
                mags += [0.0]
                errs += [0.0]
                
        return bals, mags, errs
        
    def pygetdgprimephys(fs, ph, ionic):
        cc.p_h = Q_(ph)  # set pH
        cc.ionic_strength = Q_(ionic)  # set I
        bals = []
        mags = []
        errs = []
        for f in fs:
            try:
                rxn = cc.parse_reaction_formula(f)
                isbal = rxn.is_balanced()
                bals += [isbal]
                mags += [cc.physiological_dg_prime(rxn).value.magnitude]
                errs += [cc.physiological_dg_prime(rxn).error.magnitude]
            except:
                bals += [False]
                mags += [0.0]
                errs += [0.0]
                
        return bals, mags, errs
    """
end

end # module
