"""
$(TYPEDEF)

Information about isozyme composition and activity.

# Fields
$(TYPEDFIELDS)
"""
mutable struct Isozyme
    gene_product_count::Dict{String,Int}
    kcat_forward::Float64
    kcat_reverse::Float64
end
