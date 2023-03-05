
"""
$(TYPEDEF)

A helper type that describes the contents of [`SimplifiedEnzymeConstrainedModel`](@ref)s.

# Fields
$(TYPEDFIELDS)
"""
struct SimplifiedEnzymeConstrainedColumn
    reaction_idx::Int # number of the corresponding reaction in the inner model
    direction::Int # 0 if "as is" and unique, -1 if reverse-only part, 1 if forward-only part
    lb::Float64 # must be 0 if the reaction is unidirectional (if direction!=0)
    ub::Float64
    capacity_contribution::Float64 # must be 0 for bidirectional reactions (if direction==0)
    capacity_bound_idxs::Vector{Int64} # index of associated bound(s)
end

"""
$(TYPEDEF)

A helper type for describing the contents of [`EnzymeConstrainedModel`](@ref)s.

# Fields
$(TYPEDFIELDS)
"""
struct EnzymeConstrainedReactionColumn
    reaction_idx::Int
    isozyme_idx::Int
    direction::Int
    reaction_coupling_row::Int
    lb::Float64
    ub::Float64
    gene_product_coupling::Vector{Tuple{Int,Float64}}
end

"""
$(TYPEDEF)

A helper struct that contains the gene product capacity terms organized by
the grouping type, e.g. metabolic or membrane groups etc.

# Fields
$(TYPEDFIELDS)
"""
struct EnzymeConstrainedCapacity
    group_id::String
    gene_product_idxs::Vector{Int}
    gene_product_molar_masses::Vector{Float64}
    group_upper_bound::Float64
end
