"""
$(TYPEDEF)

This is the abstract type used to map optimizer output to semantically
meaningful variables.
"""
abstract type AbstractSemanticVariables end

"""
$(TYPEDEF)

This is the concrete type that should be used to map optimizer output to 
reaction fluxes. The typical units are mmol/gDW/h.
"""
struct ReactionFluxes <: AbstractSemanticVariables end

"""
$(TYPEDEF)

This is the concrete type that should be used to map optimizer output to 
enzyme abundances. The typical units are g/gDW.
"""
struct EnzymeAbundances <: AbstractSemanticVariables end

"""
$(TYPEDEF)

This is the concrete type that should be used to map optimizer output to 
metabolite concentrations. The typical units are M.
"""
struct MetaboliteConcentrations <: AbstractSemanticVariables end
