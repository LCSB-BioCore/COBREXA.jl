function _pad_spaces(str, maxlen)
   join(repeat([" "], inner=maxlen - length(str) + 1)) 
end

function Base.show(io::IO, ::MIME"text/plain", flux_res::FluxSummary)
    longest_biomass_len = maximum([length(k) for k in keys(flux_res.biomass_fluxes)])
    println(io, "Biomass:")
    for (k, v) in flux_res.biomass_fluxes
        println(io, "  ", k, ":", _pad_spaces(k, longest_biomass_len), round(v, digits=4))
    end

    longest_import_len = maximum([length(k) for k in keys(flux_res.import_fluxes)])
    println(io, "Import:")
    for (k, v) in flux_res.import_fluxes
        println(io, "  ",  k, ":", _pad_spaces(k, longest_import_len), round(v, digits=4))
    end

    longest_export_len = maximum([length(k) for k in keys(flux_res.import_fluxes)])
    println(io, "Export:")
    for (k, v) in flux_res.export_fluxes
        println(io, "  ", k, ":", _pad_spaces(k, longest_export_len), round(v, digits=4))
    end

    if !isempty(flux_res.unbounded_fluxes)
        longest_unbounded_len = maximum([length(k) for k in keys(flux_res.unbounded_fluxes)])
        println(io, "Unbounded:")
        for (k, v) in flux_res.unbounded_fluxes
            println(io, "  ", k, ":", _pad_spaces(k, longest_unbounded_len), round(v, digits=4))
        end
    end
    return nothing
end
