"""
$(TYPEDEF)

Wrapper around the models loaded in dictionaries from the MATLAB representation.

# Fields
$(TYPEDFIELDS)
"""
struct MATModel <: MetabolicModel
    mat::Dict{String,Any}
end
