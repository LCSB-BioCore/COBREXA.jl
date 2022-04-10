"""
    _build_irreversible_stoichiometric_matrix(model::StandardModel)

Return a stoichiometric matrix where all reactions are forward only i.e. only
positive fluxes are allowed. To accomplish this for models with isozymes,
so-called arm reactions are included. Note, reactions that are irreversible 
in the original model will be irreversible in this model. E.g., if a reaction
is forward only in the original model, then there will be no reverse component 
for this reaction in the irreversible stoichiometric matrix.
"""
function _build_irreversible_stoichiometric_matrix(
    model::StandardModel,
    rid_isozymes = Dict{String,Vector{Isozyme}}(),
)
    # components used to build stoichiometric matrix
    S_components = ( #TODO add size hints if possible
        row_idxs = Vector{Int}(),
        col_idxs = Vector{Int}(),
        coeffs = Vector{Float64}(),
        lbs = Vector{Float64}(),
        ubs = Vector{Float64}(),
    )

    # establish the ordering in a named tuple
    idxs = ( #: pseudo metabolites and reactions are added to model
        met_idxs = Dict{String,Int}(),
        rxn_idxs = Dict{String,Int}(),
        max_rxn_idx = [1], #TODO maybe fix, this is a dodgy way of adding a counter to a named tuple
        max_met_idx = [1], #TODO maybe fix, this is a dodgy way of adding a counter to a named tuple
        pseudo_met_idx = [1], #TODO maybe fix, this is a dodgy way of adding a counter to a named tuple
    )
    #TODO for the counter thing, basically I wanted e.g. max_rxn_idx = 1 and then update it,
    #TODO but named tuples are immutable... :(

    # fill the matrix entries
    #: blocked treated as reversible because unclear what direction the reaction would go
    for rid in reactions(model)
        if haskey(rid_isozymes, rid) && length(rid_isozymes[rid]) > 1
            if is_reaction_unidirectional(model, rid)
                dir = is_reaction_forward_only(model, rid) ? "§FOR" : "§REV"
                _add_isozyme_to_irrev_stoich_mat(
                    model,
                    rid_isozymes,
                    rid,
                    idxs,
                    S_components,
                    dir,
                )
            elseif is_reaction_reversible(model, rid) || is_reaction_blocked(model, rid)
                _add_isozyme_to_irrev_stoich_mat(
                    model,
                    rid_isozymes,
                    rid,
                    idxs,
                    S_components,
                    "§FOR",
                )
                _add_isozyme_to_irrev_stoich_mat(
                    model,
                    rid_isozymes,
                    rid,
                    idxs,
                    S_components,
                    "§REV",
                )
            else
                @warn "Unhandled bound type for $rid"
            end
        else # no grr or single enzyme only
            if is_reaction_unidirectional(model, rid)
                dir = is_reaction_forward_only(model, rid) ? "§FOR" : "§REV"
                _add_enzyme_to_irrev_stoich_mat(model, rid, idxs, S_components, dir)
            elseif is_reaction_reversible(model, rid) || is_reaction_blocked(model, rid)
                _add_enzyme_to_irrev_stoich_mat(model, rid, idxs, S_components, "§FOR")
                _add_enzyme_to_irrev_stoich_mat(model, rid, idxs, S_components, "§REV")
            else
                @warn "Unhandled bound type for $rid"
            end
        end
    end

    S = sparse(
        S_components.row_idxs,
        S_components.col_idxs,
        S_components.coeffs,
        length(idxs.met_idxs),
        length(idxs.rxn_idxs),
    )

    return S, S_components.lbs, S_components.ubs, idxs.rxn_idxs, idxs.met_idxs
end

"""
    _add_enzyme_to_irrev_stoich_mat(model, rid, idxs, S_components, dir)

Add entries to the components that will be used to build the stoichiometric
matrix. Simple variant that does not deal with isozymes and arm reactions.
"""
function _add_enzyme_to_irrev_stoich_mat(model::StandardModel, rid, idxs, S_components, dir)
    idxs.rxn_idxs[rid*dir] = idxs.max_rxn_idx[1]
    idxs.max_rxn_idx[1] += 1
    fix_sign = dir == "§FOR" ? 1 : -1 # change direction of reaction
    for (mid, coeff) in reaction_stoichiometry(model, rid)
        if !haskey(idxs.met_idxs, mid)
            idxs.met_idxs[mid] = idxs.max_met_idx[1]
            idxs.max_met_idx[1] += 1
        end
        push!(S_components.row_idxs, idxs.met_idxs[mid])
        push!(S_components.col_idxs, idxs.rxn_idxs[rid*dir])
        push!(S_components.coeffs, fix_sign * coeff)
    end
    lb, ub = abs.(reaction_bounds(model, rid)) # assumes lb < ub
    if dir == "§FOR"
        is_reaction_reversible(model, rid) ? push!(S_components.lbs, 0) :
        push!(S_components.lbs, lb)
        push!(S_components.ubs, ub)
    else
        is_reaction_reversible(model, rid) ? push!(S_components.lbs, 0) :
        push!(S_components.lbs, ub)
        push!(S_components.ubs, lb)
    end
end

"""
    _add_isozyme_to_irrev_stoich_mat(
        model::StandardModel,
        rid,
        idxs,
        S_components,
        dir,
    )

Add entries to the components that will be used to build the stoichiometric matrix.
Complex variant that deals with isozymes and arm reactions.
"""
function _add_isozyme_to_irrev_stoich_mat(
    model::StandardModel,
    rid_isoyzmes,
    rid,
    idxs,
    S_components,
    dir,
)
    # add pseudo metabolite
    pm = "§PM$(idxs.pseudo_met_idx[1])"
    idxs.pseudo_met_idx[1] += 1
    idxs.met_idxs[pm] = idxs.max_met_idx[1]
    idxs.max_met_idx[1] += 1
    # find half reactions to get arm reaction
    lhs = []
    rhs = []
    for (mid, coeff) in reaction_stoichiometry(model, rid)
        if !haskey(idxs.met_idxs, mid)
            idxs.met_idxs[mid] = idxs.max_met_idx[1]
            idxs.max_met_idx[1] += 1
        end
        if coeff <= 0
            push!(lhs, (mid, coeff))
        else
            push!(rhs, (mid, coeff))
        end
    end
    product_half_reaction = dir == "§FOR" ? rhs : lhs
    reagent_half_reaction = dir == "§FOR" ? lhs : rhs
    # add arm reaction
    fix_sign = dir == "§FOR" ? 1 : -1 # change direction of reaction
    pr = rid * "§ARM" * dir
    idxs.rxn_idxs[pr] = idxs.max_rxn_idx[1] #! this needs to get added first because of blocked possibility
    idxs.max_rxn_idx[1] += 1
    push!(S_components.row_idxs, idxs.met_idxs[pm])
    push!(S_components.col_idxs, idxs.rxn_idxs[pr])
    push!(S_components.coeffs, 1)
    for (mid, coeff) in reagent_half_reaction
        push!(S_components.row_idxs, idxs.met_idxs[mid])
        push!(S_components.col_idxs, idxs.rxn_idxs[pr])
        push!(S_components.coeffs, fix_sign * coeff)
    end
    # add bounds for ARM reaction that corresponds to original model's bounds
    lb, ub = abs.(reaction_bounds(model, rid)) # assumes lb < ub
    if dir == "§FOR"
        is_reaction_reversible(model, rid) ? push!(S_components.lbs, 0) :
        push!(S_components.lbs, lb)
        push!(S_components.ubs, ub)
    else
        is_reaction_reversible(model, rid) ? push!(S_components.lbs, 0) :
        push!(S_components.lbs, ub)
        push!(S_components.ubs, lb)
    end
    # add isozyme reactions
    for (i, _) in enumerate(rid_isoyzmes[rid])
        iso_rid = rid * "§ISO$i" * dir
        idxs.rxn_idxs[iso_rid] = idxs.max_rxn_idx[1]
        idxs.max_rxn_idx[1] += 1
        push!(S_components.row_idxs, idxs.met_idxs[pm])
        push!(S_components.col_idxs, idxs.rxn_idxs[iso_rid])
        push!(S_components.coeffs, -1)
        for (mid, coeff) in product_half_reaction
            push!(S_components.row_idxs, idxs.met_idxs[mid])
            push!(S_components.col_idxs, idxs.rxn_idxs[iso_rid])
            push!(S_components.coeffs, fix_sign * coeff)
        end
        # add bounds
        push!(S_components.lbs, 0)
        if is_reaction_blocked(model, rid)
            push!(S_components.ubs, 0)
        else
            push!(S_components.ubs, 1000) # arbitrary upper bound
        end
    end
end

"""
    _order_id_to_idx_dict(id_to_idx_dict)

Return the keys of `id_to_idx_dict` sorted by the values, which
are taken to be the indices. This is a helper function for
[`reactions`](@ref) and [`metabolites`](@ref).
"""
function _order_id_to_idx_dict(dmap)
    ks = collect(keys(dmap))
    vs = collect(values(dmap))
    return ks[sortperm(vs)]
end
