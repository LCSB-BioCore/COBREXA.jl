function Base.show(io::IO, ::MIME"text/plain", flux_res::FluxSummary)
    c = Base.text_colors # of course I added colors :)
    # println(c[:bold]*c[:light_black], "Biomass:", c[:normal]) # display all the biomass fluxes
    # for (k, v) in zip(bmasses, bmass_fluxes)
    #     println("\t", c[:light_magenta], k, ": ", c[:normal], round(v, digits=round_digits), c[:normal])
    # end

    # println(c[:bold]*c[:light_black],"Import:", c[:normal])
    # for (k, v) in zip(ex_rxns[import_fluxes], ex_fluxes[import_fluxes])
    #     println("\t", c[:light_magenta], k, ": ", c[:normal], round(v, digits=round_digits))
    # end

    # println(c[:bold]*c[:light_black],"Export:", c[:normal])
    # for (k, v) in zip(ex_rxns[export_fluxes], ex_fluxes[export_fluxes])
    #     println("\t", c[:light_magenta], k, ": ", c[:normal], round(v, digits=round_digits))
    # end

    # println(c[:bold]*c[:light_black],"Unbounded:", c[:normal])
    # for (k, v) in zip([ex_rxns[lower_unbounded]; ex_rxns[upper_unbounded]], [ex_fluxes[lower_unbounded]; ex_fluxes[upper_unbounded]])
    #     println("\t", c[:light_magenta], k, ": ", c[:normal], round(v, digits=round_digits))
    # end

end

function Base.show(io::IO, ::MIME"text/html", flux_res::FluxSummary)

end
