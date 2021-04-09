function _get_boundary_reaction_ids(model::StandardModel)::Array{String,1}
    return [i.id for i in model.reactions if length(i.metabolites) == 1]
end

function add_cycle_free(fluxes::Dict{String,Float64})
    return (model, opt_model) -> begin
        v = opt_model[:x]
        old_objective = objective_function(opt_model)
        boundary_ids = _get_boundary_reaction_ids(model)
        min_objectives = zeros(Int64, 1, length(v))
        for (i, reaction) in enumerate(model.reactions)
        id = reaction.id
        if v[i] in old_objective.terms.keys
            continue
        end
        flux = fluxes[id]
        if id in boundary_ids
            lb = flux
            ub = flux
        elseif flux >= 0
                lb = max(0, reaction.lb)
                ub = max(flux, reaction.ub)
                min_objectives[i] =  1
            else
                lb = min(flux, reaction.lb)
                ub = min(0, reaction.ub)            
                min_objectives[i] = -1
        end
        set_bound(i, opt_model; ub=ub, lb=lb)
    end    
        @objective(opt_model, Min, dot(min_objectives, v))
    end
end


"""
    add_loopless(args...)::Union{COBREXA.JuMP.Model,Nothing}

Convert an existing FBA solution to a loopless one. Removes as many loops as possible using the method from [1]

References:
[1] CycleFreeFlux: efficient removal of thermodynamically infeasible
    loops from flux distributions. Desouki AA, Jarre F, Gelius-Dietrich
    G, Lercher MJ. Bioinformatics. 2015 Jul 1;31(13):2159-65. doi:
    10.1093/bioinformatics/btv096.
"""
function add_loopless(model::M, opt_model, fluxes::Dict{String,Float64})::Union{COBREXA.JuMP.Model,Nothing} where {M <: MetabolicModel}
    constrain_objective_value(1)(model, opt_model)
    add_cycle_free(fluxes)(model, opt_model)
    COBREXA.JuMP.optimize!(opt_model)    
    termination_status(opt_model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED] || return nothing
    return opt_model
end

"""
    add_loopless_vec(args...)::Union{Array{Float64,1},Nothing}

Convert an existing FBA solution to a loopless one. Removes as many loops as possible using the method from [1]

References:
[1] CycleFreeFlux: efficient removal of thermodynamically infeasible
    loops from flux distributions. Desouki AA, Jarre F, Gelius-Dietrich
    G, Lercher MJ. Bioinformatics. 2015 Jul 1;31(13):2159-65. doi:
    10.1093/bioinformatics/btv096.
"""
function add_loopless_vec(model::M, opt_model, fluxes::Dict{String,Float64})::Union{Array{Float64,1},Nothing} where {M <: MetabolicModel}
    opt_model = add_loopless(model, opt_model, fluxes)
    termination_status(opt_model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED] || return nothing
    return value.(opt_model[:x])
end

"""
    add_loopless_dict(args...)::Union{Dict{String,Float64},Nothing}

Convert an existing FBA solution to a loopless one. Removes as many loops as possible using the method from [1]

References:
[1] CycleFreeFlux: efficient removal of thermodynamically infeasible
    loops from flux distributions. Desouki AA, Jarre F, Gelius-Dietrich
    G, Lercher MJ. Bioinformatics. 2015 Jul 1;31(13):2159-65. doi:
    10.1093/bioinformatics/btv096.
"""
function add_loopless_dict(model::M, opt_model, fluxes::Dict{String,Float64})::Union{Dict{String,Float64},Nothing} where {M <: MetabolicModel}
    v = add_loopless_vec(model, opt_model, fluxes)
    isnothing(v) && return nothing
    return Dict(zip(reactions(model), v))
end

"""
    loopless_flux_balance_analysis(args...)::Union{COBREXA.JuMP.Model,Nothing}

Run loopless FBA. Removes as many loops as possible using the method from [1]

References:
[1] CycleFreeFlux: efficient removal of thermodynamically infeasible
    loops from flux distributions. Desouki AA, Jarre F, Gelius-Dietrich
    G, Lercher MJ. Bioinformatics. 2015 Jul 1;31(13):2159-65. doi:
    10.1093/bioinformatics/btv096.
"""
function loopless_flux_balance_analysis(
    model::M,
    optimizer;
    modifications=[(model, opt_model) -> nothing],

)::Union{COBREXA.JuMP.Model,Nothing} where {M <: MetabolicModel}
    opt_model = flux_balance_analysis(model, optimizer,  modifications=modifications)
    termination_status(opt_model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED] || return nothing
    fluxes = Dict(zip(reactions(model), value.(opt_model[:x])))
    return add_loopless(
        model,
        opt_model,
        fluxes,
    )
end

"""
    loopless_flux_balance_analysis_vec(args...)::Union{Array{Float64,1},Nothing}

Run loopless FBA. Removes as many loops as possible using the method from [1]

References:
[1] CycleFreeFlux: efficient removal of thermodynamically infeasible
    loops from flux distributions. Desouki AA, Jarre F, Gelius-Dietrich
    G, Lercher MJ. Bioinformatics. 2015 Jul 1;31(13):2159-65. doi:
    10.1093/bioinformatics/btv096.
"""
function loopless_flux_balance_analysis_vec(
    model::M,
    optimizer;
    modifications=[(model, opt_model) -> nothing],

)::Union{Array{Float64,1},Nothing} where {M <: MetabolicModel}
    opt_model = loopless_flux_balance_analysis(model, optimizer, modifications=modifications)
    termination_status(opt_model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED] || return nothing
    return value.(opt_model[:x])
end

"""
    loopless_flux_balance_analysis_dict(args...)::Union{Dict{String,Float64},Nothing}

Run loopless FBA. Removes as many loops as possible using the method from [1]

References:
[1] CycleFreeFlux: efficient removal of thermodynamically infeasible
    loops from flux distributions. Desouki AA, Jarre F, Gelius-Dietrich
    G, Lercher MJ. Bioinformatics. 2015 Jul 1;31(13):2159-65. doi:
    10.1093/bioinformatics/btv096.
"""
function loopless_flux_balance_analysis_dict(
    model::M,
    optimizer;
    modifications=[(model, opt_model) -> nothing],

)::Union{Dict{String,Float64},Nothing} where {M <: MetabolicModel}
    v = loopless_flux_balance_analysis_vec(model, optimizer, modifications=modifications)
    isnothing(v) && return nothing
    return Dict(zip(reactions(model), v))
end
