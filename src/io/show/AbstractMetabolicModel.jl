"""
$(TYPEDSIGNATURES)

Pretty printing of everything metabolic-modelish.
"""
function Base.show(io::Base.IO, ::MIME"text/plain", m::AbstractMetabolicModel)
    _pretty_print_keyvals(
        io,
        "",
        "$(typeof(m)){$(n_reactions(m)) reactions, $(n_metabolites(m)) metabolites}",
    )
end
