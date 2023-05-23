function Base.show(io::Base.IO, ::MIME"text/plain", res::ModelWithResult{T}) where {T}
    println(io, "ModelWithResult{$T}($(typeof(res.model)), ...)")
end