"""
    atom_exchange(flux_dict::Dict{String, Float64}, model::StandardModel)

Return a dictionary mapping the flux of atoms across the boundary of the model 
given `flux_dict` (the solution of a constraint based analysis) of reactions in `model`.
"""
function atom_exchange(flux_dict::Dict{String,Float64}, model::StandardModel)
    atom_flux = Dict{String,Float64}()
    for (rxn_id, flux) in flux_dict
        if is_boundary(model.reactions[rxn_id])
            for (met, stoich) in model.reactions[rxn_id].metabolites
                adict = get_atoms(model.metabolites[met])
                for (atom, stoich) in adict
                    atom_flux[atom] = get(atom_flux, atom, 0.0) + flux * stoich
                end
            end
        end
    end
    return atom_flux
end

"""
    get_exchanges(flux_dict::Dict{String, Float64}; top_n=Inf, ignorebound=_constants.default_reaction_bound, verbose=true)

Display the `top_n` producing and consuming exchange fluxes.
If `top_n` is not specified (by an integer), then all are displayed.
Ignores infinite (problem upper/lower bound) fluxes (set with ignorebound).
When `verbose` is false, the output is not printed out.
Return these reactions (id => ) in two dictionaries: `consuming`, `producing`
"""
function exchange_reactions(
    flux_dict::Dict{String,Float64},
    model::StandardModel;
    top_n = Inf,
    ignorebound = _constants.default_reaction_bound,
    verbose = true,
)
    consuming = Dict{String,Float64}()
    producing = Dict{String,Float64}()

    for (k, v) in flux_dict
        if is_boundary(model.reactions[k])
            if v < 0 # consuming
                consuming[k] = v
            elseif v > 0 # producing
                producing[k] = v
            else # no flux
                continue
            end
        end
    end

    if verbose
        # Do consuming
        ks = collect(keys(consuming))
        vs = [consuming[k] for k in ks]
        inds = sortperm(vs)
        n_max = length(ks)
        println("Consuming fluxes: ")
        ii = 0 # counter
        for i in inds
            if vs[i] > -ignorebound
                println(ks[i], " = ", round(vs[i], digits = 6))
                ii += 1
            end
            if ii > top_n
                break
            end
        end
        # Do producing
        ks = collect(keys(producing))
        vs = [producing[k] for k in ks]
        inds = sortperm(vs, rev = true)
        n_max = length(ks)
        println("Producing fluxes: ")
        ii = 0 # counter
        for i in inds
            if vs[i] < ignorebound
                println(ks[i], " = ", round(vs[i], digits = 6))
                ii += 1
            end
            if ii > top_n
                break
            end
        end
    end

    return consuming, producing
end

"""
    metabolite_fluxes(flux_dict::Dict{String, Float64}, model::StandardModel)

Return two dictionaries of metabolite `id`s mapped to reactions that consume or 
produce them given the flux distribution supplied in `fluxdict`.
"""
function metabolite_fluxes(flux_dict::Dict{String,Float64}, model::StandardModel)
    S = stoichiometry(model)
    met_flux = Dict{String,Float64}()
    rxnids = reactions(model)
    metids = metabolites(model)

    producing = Dict{String,Dict{String,Float64}}()
    consuming = Dict{String,Dict{String,Float64}}()
    for (row, metid) in enumerate(metids)
        for (col, rxnid) in enumerate(rxnids)
            mf = flux_dict[rxnid] * S[row, col]
            # ignore zero flux
            if mf < -_constants.tolerance # consuming rxn
                if haskey(consuming, metid)
                    consuming[metid][rxnid] = mf
                else
                    consuming[metid] = Dict(rxnid => mf)
                end
            elseif mf > _constants.tolerance
                if haskey(producing, metid)
                    producing[metid][rxnid] = mf
                else
                    producing[metid] = Dict(rxnid => mf)
                end
            end
        end
    end
    return consuming, producing
end

"""
    set_bound(index, optimization_model;
        ub=_constants.default_reaction_rate,
        lb=-_constants.default_reaction_rate)

Helper function to set the bounds of variables.
The JuMP `set_normalized_rhs` function is a little confusing, 
so this function simplifies setting constraints. In short, JuMP
uses a normalized right hand side representation of constraints, 
which means that lower bounds have their sign flipped. This function
does this for you, so you don't have to remember to do this whenever you
change the constraints. 

Just supply the constraint `index` and the JuMP model (`opt_model`) that 
will be solved, and the variable's bounds will be set to `ub` and `lb`.
"""
function set_bound(
    vind,
    opt_model;
    ub = _constants.default_reaction_rate,
    lb = -_constants.default_reaction_rate,
)
    set_normalized_rhs(opt_model[:lbs][vind], -lb)
    set_normalized_rhs(opt_model[:ubs][vind], ub)
end

"""
    get_bound_vectors(opt_model)

Returns vectors of the lower and upper bounds of `opt_model` constraints, where
`opt_model` is a JuMP model constructed by e.g.
[`make_optimization_model`](@ref) or [`flux_balance_analysis`](@ref).


"""
function get_bound_vectors(opt_model)
    lbconref = opt_model[:lbs]
    ubconref = opt_model[:ubs]
    lbs = zeros(length(lbconref))
    for i in eachindex(lbs)
        lbs[i] = -normalized_rhs(lbconref[i])
    end
    ubs = [normalized_rhs(ubconref[i]) for i in eachindex(ubconref)]

    return lbs, ubs
end
