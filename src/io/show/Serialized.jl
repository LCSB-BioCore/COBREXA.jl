
"""
$(TYPEDSIGNATURES)

Show the [`Serialized`](@ref) model without unnecessarily loading it.
"""
function Base.show(io::IO, ::MIME"text/plain", m::Serialized{M}) where {M}
    print(
        io,
        "Serialized{$M} saved in \"$(m.filename)\" ($(isnothing(m.m) ? "not loaded" : "loaded"))",
    )
end
