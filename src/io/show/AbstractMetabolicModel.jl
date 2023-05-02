"""
$(TYPEDSIGNATURES)

Pretty printing of everything metabolic-modelish.
"""
function Base.show(io::Base.IO, ::MIME"text/plain", m::AbstractMetabolicModel)
    print(
        io,
        "$(typeof(m))(#= $(reaction_count(m)) reactions, $(n_metabolites(m)) metabolites =#)",
    )
end
