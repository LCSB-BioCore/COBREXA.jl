"""
$(TYPEDSIGNATURES)

Pretty printing of a [`AbstractResult`](@ref) and more.
"""
function Base.show(io::Base.IO, ::MIME"text/plain", m::AbstractResult)
    print(
        io,
        """
        Result with fields `model` and `opt_model`.
        """,
    )
end

# $(typeof(m))(#= $(n_reactions(m)) reactions, $(n_metabolites(m)) metabolites =#)\n 
