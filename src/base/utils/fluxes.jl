"""
    metabolite_fluxes(model::MetabolicModel, flux_dict::Dict{String, Float64})

Return two dictionaries of metabolite `id`s mapped to reactions that consume or
produce them, given the flux distribution supplied in `flux_dict`.
"""
function metabolite_fluxes(model::MetabolicModel, flux_dict::Dict{String,Float64})
    S = stoichiometry(model)
    rids = reactions(model)
    mids = metabolites(model)

    producing = Dict{String,Dict{String,Float64}}()
    consuming = Dict{String,Dict{String,Float64}}()
    for (row, mid) in enumerate(mids)
        for (col, rid) in enumerate(rids)
            mf = flux_dict[rid] * S[row, col]
            # ignore zero flux
            if mf < -_constants.tolerance # consuming rxn
                if haskey(consuming, mid)
                    consuming[mid][rid] = mf
                else
                    consuming[mid] = Dict(rid => mf)
                end
            elseif mf > _constants.tolerance
                if haskey(producing, mid)
                    producing[mid][rid] = mf
                else
                    producing[mid] = Dict(rid => mf)
                end
            end
        end
    end
    return consuming, producing
end
