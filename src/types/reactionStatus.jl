"""
Used for concise reporting of modeling results.
"""
mutable struct ReactionStatus
    alreadyPresent::Bool
    index::Int
    info::String
end
