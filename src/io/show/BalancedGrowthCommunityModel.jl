function Base.show(io::Base.IO, ::MIME"text/plain", cm::CommunityMember)
    println(io, "A $(typeof(cm.model)) CommunityMember with:")
    println(io, "  id $(cm.id),")
    println(io, "  abundance $(cm.abundance),")
    println(
        io,
        "  underlying model of $(n_reactions(cm.model)) reactions and $(n_metabolites(cm.model)) metabolites,",
    )
    println(io, "  biomass metabolite $(cm.biomass_metabolite_id),")
    println(
        io,
        "  and $(length(cm.exchange_reaction_ids)) exchange reactions that will be connected to the environment.",
    )
end

function Base.show(io::Base.IO, ::MIME"text/plain", cm::BalancedGrowthCommunityModel)
    println(io, "A balanced growth community model with:")
    println(io, "  $(length(cm.members)) underlying models,")
    println(io, "  objective $(cm.objective_id),")
    if isempty(cm.env_met_flux_bounds)
        println(io, "  and no constraints on environmental metabolite fluxes.")
    else
        println(io, "  and constraints on the following environmental metabolite fluxes:")
        for (k, v) in cm.env_met_flux_bounds
            println(io, "    $(first(v)) ≤ $(k) ≤ $(last(v))")
        end
    end
end
