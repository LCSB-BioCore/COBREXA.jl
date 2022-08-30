"""
$(TYPEDEF)

Used for concise reporting of modeling results.

# Fields
$(TYPEDFIELDS)
"""
mutable struct ReactionStatus
    already_present::Bool
    index::Int
    info::String
end
