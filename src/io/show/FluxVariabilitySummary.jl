function Base.show(io::Base.IO, ::MIME"text/plain", flux_res::FluxVariabilitySummary)

    longest_biomass_len =
        maximum(length(k) for k in keys(flux_res.biomass_fluxes); init = 0)
    longest_exchange_len =
        maximum(length(k) for k in keys(flux_res.exchange_fluxes); init = 0)
    word_pad_len = max(longest_biomass_len, longest_exchange_len)

    longest_biomass_flux_len = maximum(
        length(string(round(v[1], digits = 4))) for v in values(flux_res.biomass_fluxes);
        init = 0,
    )
    longest_exchange_flux_len = max(
        maximum(
            length(string(round(v[1], digits = 4))) for
            v in values(flux_res.exchange_fluxes);
            init = 0,
        ),
        length("Lower bound  "),
    )
    number_pad_len = max(longest_biomass_flux_len, longest_exchange_flux_len)

    println(
        io,
        "Biomass",
        _pad_spaces(length("Biomass"), word_pad_len + 3),
        "Lower bound  ",
        _pad_spaces("Lower bound  ", number_pad_len),
        "Upper bound",
    )
    for (k, v) in flux_res.biomass_fluxes
        lb = isnothing(v[1]) ? " " : string(round(v[1], digits = 4))
        ub = isnothing(v[2]) ? " " : string(round(v[1], digits = 4))
        println(
            io,
            "  ",
            k,
            ":",
            _pad_spaces(k, word_pad_len),
            lb,
            _pad_spaces(lb, number_pad_len),
            ub,
        )
    end

    println(io, "Exchange")
    ex_ids = collect(keys(flux_res.exchange_fluxes))
    vs = [
        abs(flux_res.exchange_fluxes[k][1]) + abs(flux_res.exchange_fluxes[k][2]) for
        k in ex_ids
    ]
    idxs = sortperm(vs, rev = true)
    for k in ex_ids[idxs]
        v = flux_res.exchange_fluxes[k]
        lb = isnothing(v[1]) ? " " : string(round(v[1], digits = 4))
        ub = isnothing(v[2]) ? " " : string(round(v[2], digits = 4))
        println(
            io,
            "  ",
            k,
            ":",
            _pad_spaces(k, word_pad_len),
            lb,
            _pad_spaces(lb, number_pad_len),
            ub,
        )
    end
end
