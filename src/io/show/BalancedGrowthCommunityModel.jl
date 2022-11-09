function Base.show(io::Base.IO, ::MIME"text/plain", cm::CommunityMember)
    println(io, "A $(typeof(cm.model)) community member with $(n_reactions(cm.model)) reactions, $(n_metabolites(cm.model)) metabolites, and abundance $(cm.abundance*100)%.")
end

function Base.show(io::Base.IO, ::MIME"text/plain", cm::BalancedGrowthCommunityModel)
    println(io, "A balanced growth community model comprised of $(length(cm.members)) underlying models.")
end
