"""
Tests if two `StandardModel`'s represent the same model core model internally.
"""
function Base.isequal(model1::StandardModel, model2::StandardModel)
    # test if blank model is given - automatic fail
    (isempty(model1.reactions) || isempty(model2.reactions)) ? (return false) : nothing

    # test same rxn and met ids
    rxns1 = reactions(model1)
    rxns2 = reactions(model2)
    mets1 = metabolites(model1)
    mets2 = metabolites(model2)
    rxns_same =
        (isempty(setdiff(rxns1, rxns2)) && isempty(setdiff(rxns2, rxns1))) ? true : false
    mets_same =
        (isempty(setdiff(mets1, mets2)) && isempty(setdiff(mets2, mets1))) ? true : false

    if !rxns_same || !mets_same # if the ids are different stop testing
        return false
    end

    for rxn_id in rxns1 # since IDs are the same only use the ids of the first model
        # check stoichiometry
        rmets1 = model1.reactions[rxn_id].metabolites
        rmets2 = model2.reactions[rxn_id].metabolites
        if length(rmets1) != length(rmets2)
            return false
        end
        for (km, vm) in rmets1
            if !(km in keys(rmets2))
                return false
            elseif rmets2[km] != vm
                return false
            end
        end

        # check bounds
        if model1.reactions[rxn_id].lower_bound != model2.reactions[rxn_id].lower_bound
            return false
        end
        if model1.reactions[rxn_id].upper_bound != model2.reactions[rxn_id].upper_bound
            return false
        end
    end

    return true
end

"""
    atom_exchange(flux_dict::Dict{String, Float64}, model::StandardModel)

Return a dictionary mapping the flux of atoms across the boundary of the model given `flux_dict` of reactions in `model`.
Here `flux_dict` is a mapping of reaction `id`s to fluxes, e.g. from FBA.
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
    get_exchanges(rxndict::Dict{String, Float64}; top_n=Inf, ignorebound=_constants.default_reaction_bound, verbose=true)

Display the top_n producing and consuming exchange fluxes.
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
            if v[i] > -ignorebound
                println(ks[i], " = ", round(v[i], digits = 6))
                ii += 1
            end
            if ii > top_n
                break
            end
        end
        # Do producing
        ks = collect(keys(producing))
        vs = [producing[k] for k in ks]
        inds = sortperm(vs)
        n_max = length(ks)
        println("Producing fluxes: ")
        ii = 0 # counter
        for i in inds
            if v[i] < ignorebound
                println(ks[i], " = ", round(v[i], digits = 6))
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
    metabolite_fluxes(fluxdict::Dict{String, Float64}, model::StandardModel)

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
    if lb <= 0
        set_normalized_rhs(opt_model[:lbs][vind], abs(lb))
    else
        set_normalized_rhs(opt_model[:lbs][vind], -abs(lb))
    end
    set_normalized_rhs(opt_model[:ubs][vind], ub)
end

"""
    get_bound_vectors(opt_model)

Returns vectors of the lower and upper bounds of `opt_model` constraints, where 
`opt_model` is a JuMP model constructed by e.g. `make_optimization_problem` or
`flux_balance_analysis`.

See also: [`make_optimization_problem`](@ref), [`flux_balance_analysis`](`ref`)
"""
function get_bound_vectors(opt_model)
    lbconref = opt_model[:lbs]
    ubconref = opt_model[:ubs]
    lbs = zeros(length(lbconref))
    for i in eachindex(lbs)
        lbval = normalized_rhs(lbconref[i])
        if lbval > 0
            lbs[i] = -abs(lbval)
        else
            lbs[i] = abs(lbval)
        end
    end
    ubs = [normalized_rhs(ubconref[i]) for i in eachindex(ubconref)]

    return lbs, ubs
end

"""
    is_mass_balanced(rxn::Reaction, model::StandardModel)

Checks if `rxn` is atom balanced. Returns a boolean for whether the reaction is balanced,
and the associated balance of atoms for convenience (useful if not balanced).

See also: [`get_atoms`](@ref), [`check_duplicate_reaction`](@ref)
"""
function is_mass_balanced(rxn::Reaction, model::StandardModel)
    atom_balances = Dict{String,Float64}() # float here because stoichiometry is not Int
    for (met, stoich) in rxn.metabolites
        atoms = get_atoms(model.metabolites[met])
        isempty(atoms) && continue # ignore blanks
        for (k, v) in atoms
            atom_balances[k] = get(atom_balances, k, 0) + v * stoich
        end
    end

    return all(sum(values(atom_balances)) == 0), atom_balances
end

"""
    knockout!(model::StandardModel, gene_id::String)

Knockout a gene in the model, which will remove all affected reactions
"""
function knockout!(model::StandardModel, gene_id::String)
    # Dev note: the three nested for loops are inefficiency. However:
    # - gene_ids (user input) will be probably only very few items
    # - model.genes[gene_id].reactions are just a few reactions (most genes don't code for a lot of reactions)
    # - reaction.grr also should only hold few items (reactions aren't coded by many different combinations of genes)
    # Let's avoid premature optimization for now and see if anyone ever has problems with this
    for reaction_id in gene_associated_reactions(model, gene_id)
        reaction = model.reactions[reaction_id]
        for (i, gene_array) in enumerate(reaction.grr)
            for gene in gene_array
                # AND inside the gene_array, so destroy as soon as one is missing
                if gene == gene_id
                    deleteat!(reaction.grr, i)
                    break
                end
            end

            # OR outside, so all have to be deleted for the reaction to be deleted
            if length(reaction.grr) == 0
                rm!(Reaction, model, reaction.id)
            end
        end
    end
    return nothing
end
"""
    knockout!(model::StandardModel, gene_ids::Vector{String})

Knockout genes in the model, which will remove all affected reactions
"""
function knockout!(model::StandardModel, gene_ids::Vector{String})
    for gene_id in gene_ids
        knockout!(model, gene_id)
    end
end
