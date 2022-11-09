"""
$(TYPEDSIGNATURES)

Pretty printing of everything metabolic-modelish.
"""
function Base.show(io::Base.IO, ::MIME"text/plain", m::AbstractMetabolicModel)
    _pretty_print_keyvals(io, "", "Metabolic model of type $(typeof(m)) with $(n_reactions(m)) reactions and $(n_metabolites(m)) metabolites.")
end
