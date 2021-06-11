function Base.show(io::IO, ::MIME"text/plain", flux_res::FluxVariabilitySummary)

    longest_biomass_len = maximum([length(k) for k in keys(flux_res.biomass_fluxes)])
    longest_exchange_len = maximum([length(k) for k in keys(flux_res.exchange_fluxes)])
    word_pad_len = max(longest_biomass_len, longest_exchange_len)

    longest_biomass_flux_len = maximum([length(string(round(v[1], digits)) for v in values(flux_res.biomass_fluxes)])
    longest_exchange_flux_len = maximum([length(string(v[1])) for v in values(flux_res.exchange_fluxes)])
    number_pad_len = max(longest_biomass_flux_len, longest_exchange_flux_len)

    println(io, "Biomass:", _pad_spaces(length("Biomass:  "), word_pad_len),"Lower bound", _pad_spaces("Lower bound", number_pad_len), "Upper bound")
    for (k, v) in flux_res.biomass_fluxes
        lb = isnothing(v[1]) ? " " : string(round(v[1], digits=4))
        ub = isnothing(v[2]) ? " " : string(round(v[1], digits=4))
        println(io, "  ", k, ":", _pad_spaces(k, longest_biomass_len), lb, _pad_spaces(lb, longest_biomass_flux_len), ub)
    end

    # println(io, "Exchange",join(repeat([" "], inner=longest_exchange_len+2-6)),"Lower bound", _pad_spaces("Lower bound", longest_exchange_flux_len), "Upper bound")
    # for (k, v) in flux_res.exchange_fluxes
    #     lb = isnothing(v[1]) ? " " : string(round(v[1], digits=4))
    #     ub = isnothing(v[2]) ? " " : string(round(v[1], digits=4))
    #     println(io, "  ", k, ":", _pad_spaces(k, longest_exchange_len), lb, _pad_spaces(lb, longest_exchange_flux_len), ub)
    # end
end
