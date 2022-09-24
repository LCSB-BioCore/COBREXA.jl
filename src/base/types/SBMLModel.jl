"""
$(TYPEDEF)

Thin wrapper around the model from SBML.jl library. Allows easy conversion from
SBML to any other model format.

# Fields
$(TYPEDFIELDS)
"""
struct SBMLModel <: MetabolicModel
    sbml::SBML.Model
end
