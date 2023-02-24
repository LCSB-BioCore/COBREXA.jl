function Base.show(io::Base.IO, ::MIME"text/plain", res::ModelWithResult)
    println(io, "A ModelWithResult composed of:")
    println(io, "model: $(typeof(res.model))")
    println(io, "result: $(typeof(res.result))")
end
