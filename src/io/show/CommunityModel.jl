function Base.show(io::Base.IO, ::MIME"text/plain", cm::CommunityMember)
    println(io, "A community member with $(variable_count(cm.model)) internal variables.")
end

function Base.show(io::Base.IO, ::MIME"text/plain", cm::CommunityModel)
    println(
        io,
        "A basic community model comprised of $(length(cm.members)) underlying models.",
    )
end
