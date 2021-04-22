"""
Pretty printing of `Gene`.
"""
function Base.show(io::IO, ::MIME"text/plain", g::Gene)
    println(io, "Gene ID: ", g.id)
    println(io, "Gene name: ", g.name)
end

"""
Pretty printing of `Vector{Gene}`.
"""
function Base.show(io::IO, ::MIME"text/plain", gs::Vector{Gene})
    println(io, "Gene set of length: ", length(gs))
end
