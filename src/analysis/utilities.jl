"""
    atom_exchange(flux_dict::Dict{String, Float64}, model::StandardModel)

Return a dictionary mapping the flux of atoms across the boundary of the model given `flux_dict` of reactions in `model`.
Here `flux_dict` is a mapping of reaction `id`s to fluxes, e.g. from FBA.
"""
function atom_exchange(flux_dict::Dict{String,Float64}, model::StandardModel)
    atom_flux = Dict{String,Float64}()
    for (rxnid, flux) in flux_dict
        if startswith(rxnid, "EX_") || startswith(rxnid, "DM_") # exchange, demand reaction
            for (met, stoich) in findfirst(model.reactions, rxnid).metabolites
                adict = get_atoms(met)
                for (atom, stoich) in adict
                    atom_flux[atom] = get(atom_flux, atom, 0.0) + flux * stoich
                end
            end
        end
    end
    return atom_flux
end

"""
    get_exchanges(rxndict::Dict{String, Float64}; top_n=8, ignorebound=1000.0, verbose=true)

Display the top_n producing and consuming exchange fluxes.
Set top_n to a large number to get all the consuming/producing fluxes.
Ignores infinite (problem upper/lower bound) fluxes (set with ignorebound).
When `verbose` is false, the output is not printed out.
Return these reactions in two dictionaries: `consuming`, `producing`
"""
function exchange_reactions(
    rxndict::Dict{String,Float64};
    top_n = 8,
    ignorebound = 1000.0,
    verbose = true,
)
    fluxes = Float64[]
    rxns = String[]
    for (k, v) in rxndict
        if startswith(k, "EX_") && abs(v) < ignorebound
            push!(rxns, k)
            push!(fluxes, v)
        end
    end
    inds_prod = sortperm(fluxes, rev = true)
    inds_cons = sortperm(fluxes)

    consuming = Dict{String,Float64}()
    producing = Dict{String,Float64}()
    verbose && println("Consuming fluxes:")
    for i = 1:min(top_n, length(rxndict))
        if rxndict[rxns[inds_cons[i]]] < -eps()
            verbose && println(
                rxns[inds_cons[i]],
                " = ",
                round(rxndict[rxns[inds_cons[i]]], digits = 4),
            )
            consuming[rxns[inds_cons[i]]] = rxndict[rxns[inds_cons[i]]]
        else
            continue
        end
    end

    verbose && println("Producing fluxes:")
    for i = 1:min(top_n, length(rxndict))
        if rxndict[rxns[inds_prod[i]]] > eps()
            verbose && println(
                rxns[inds_prod[i]],
                " = ",
                round(rxndict[rxns[inds_prod[i]]], digits = 4),
            )
            producing[rxns[inds_prod[i]]] = rxndict[rxns[inds_prod[i]]]
        else
            continue
        end
    end
    return consuming, producing
end

"""
    metabolite_fluxes(fluxdict::Dict{String, Float64}, model::StandardModel)

Return two dictionaries of metabolite `id`s mapped to reactions that consume or produce them given the flux distribution supplied in `fluxdict`.
"""
function metabolite_fluxes(fluxdict::Dict{String,Float64}, model::StandardModel)
    S = Array(stoichiometry(model))
    met_flux = Dict{String,Float64}()
    rxnids = reactions(model)
    metids = metabolites(model)

    producing = Dict{String,Dict{String,Float64}}()
    consuming = Dict{String,Dict{String,Float64}}()
    for (row, metid) in enumerate(metids)
        for (col, rxnid) in enumerate(rxnids)
            mf = fluxdict[rxnid] * S[row, col]
            # ignore zero flux
            if mf < -eps() # consuming rxn
                if haskey(consuming, metid)
                    consuming[metid][rxnid] = mf
                else
                    consuming[metid] = Dict(rxnid => mf)
                end
            elseif mf > eps()
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
    set_bound(index, optimization_model; ub=1000, lb=-1000)
Helper function to set the bounds of variables.
The JuMP `set_normalized_rhs` function is a little confusing, so this function simplifies setting
constraints. Just supply the constraint `index` and the model and they will be set to `ub` and `lb`.
"""
function set_bound(vind, opt_model; ub = 1000, lb = -1000)
    if lb <= 0
        set_normalized_rhs(opt_model[:lbs][vind], abs(lb))
    else
        set_normalized_rhs(opt_model[:lbs][vind], -abs(lb))
    end
    set_normalized_rhs(opt_model[:ubs][vind], ub)
end

"""
    modify_constraint(reaction::Reaction, lb, ub)

Modify constraints of model reaction.
"""
function modify_constraint(reaction::Reaction, lb, ub)
    (model, opt_model) -> begin
        ind = model.reactions[reaction]
        set_bound(ind, opt_model, lb = lb, ub = ub)
    end
end

"""
    modify_solver_attribute(option_key, option_val)

Modify a solver attribute. These attributes are solver specific,
refer the either JuMP or the solver you are using's documentation.
"""
function modify_solver_attribute(option_key, option_val)
    (model, opt_model) -> begin
        JuMP.set_optimizer_attribute(opt_model, option_key, option_val)
    end
end

"""
    modify_sense(objective_sense)

Modify the objective sense. 
Possible arguments are `MOI.MAX_SENSE` and `MOI.MIN_SENSE`.
Note, [`modify_objective`](@ref) sets the sense of the objective,
so it doesn't make sense to use this function AND [`modify_objective`](@ref) simultaneously. 
"""
function modify_sense(objective_sense)
    (model, opt_model) -> begin
        JuMP.set_objective_sense(opt_model, objective_sense)
    end
end

"""
    modify_objective(objective_functions::Union{Reaction, Array{Reaction, 1}}; weights=[], sense=MOI.MAX_SENSE)

Callback function to modify the objective function used in a constraint based analysis function. 
`objective_functions` can be a single reaction or an array of reactions (of type `Reaction`).
Optionally specify their `weights` and the sense of the new objective (`MOI.MAX_SENSE`, `MOI.MIN_SENSE`).
Note, this function sets the sense of the objective.
"""
function modify_objective(
    objective_functions::Union{Reaction,Array{Reaction,1}};
    weights = [],
    sense = MOI.MAX_SENSE,
)
    (model, opt_model) -> begin

        # Construct objective_indices array
        if typeof(objective_functions) == Reaction
            objective_indices = [model[objective_functions]]
        else
            objective_indices = [model[rxn] for rxn in objective_functions]
        end

        # Initialize weights
        opt_weights = zeros(length(model.reactions))

        isempty(weights) && (weights = ones(length(objective_indices))) # equal weights

        wcounter = 1
        for i in eachindex(model.reactions)
            if i in objective_indices
                opt_weights[i] = weights[wcounter]
                wcounter += 1
            end
        end

        v = opt_model[:x]
        @objective(opt_model, sense, sum(opt_weights[i] * v[i] for i in objective_indices))
    end
end

"""
    modify_solver(optimizer)

Modify the solver used to solve the model.
Typically the solver is specified as a required argument. 
This function is useful if the problem has multiple subparts that require different solvers.
See 
"""
function modify_solver(optimizer)
    (model, opt_model) -> begin
        JuMP.set_optimizer(opt_model, optimizer)
    end
end
