function _pad_spaces(slen::Int, maxlen::Int)
    " "^max(0, maxlen - slen + 1)
end

_pad_spaces(str::String, maxlen::Int) = _pad_spaces(length(str), maxlen)

function Base.show(io::IO, ::MIME"text/plain", flux_res::FluxSummary)
    longest_biomass_len =
        maximum(length(k) for k in keys(flux_res.biomass_fluxes); init = 0)
    longest_import_len = maximum(length(k) for k in keys(flux_res.import_fluxes); init = 0)
    longest_export_len = maximum(length(k) for k in keys(flux_res.export_fluxes); init = 0)

    if !isempty(flux_res.unbounded_fluxes)
        longest_unbounded_len =
            maximum([length(k) for k in keys(flux_res.unbounded_fluxes)])
        word_pad = max(
            longest_biomass_len,
            longest_export_len,
            longest_import_len,
            longest_unbounded_len,
        )
    else
        word_pad = max(longest_biomass_len, longest_export_len, longest_import_len)
    end

    println(io, "Biomass")
    for (k, v) in flux_res.biomass_fluxes
        println(io, "  ", k, ":", _pad_spaces(k, word_pad), round(v, digits = 4))
    end

    println(io, "Import")
    for (k, v) in flux_res.import_fluxes
        println(io, "  ", k, ":", _pad_spaces(k, word_pad), round(v, digits = 4))
    end

    println(io, "Export")
    for (k, v) in flux_res.export_fluxes
        println(io, "  ", k, ":", _pad_spaces(k, word_pad), round(v, digits = 4))
    end

    if !isempty(flux_res.unbounded_fluxes)
        println(io, "Unbounded")
        for (k, v) in flux_res.unbounded_fluxes
            println(io, "  ", k, ":", _pad_spaces(k, word_pad), round(v, digits = 4))
        end
    end
    return nothing
end
