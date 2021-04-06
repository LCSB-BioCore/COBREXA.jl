Base.isequal(model1::LinearModel, model2::LinearModel) =
    isequal(model1.S, model2.S) &&
    isequal(model1.b, model2.b) &&
    isequal(model1.c, model2.c) &&
    isequal(model1.xl, model2.xl) &&
    isequal(model1.xu, model2.xu) &&
    isequal(model1.rxns, model2.rxns) &&
    isequal(model1.mets, model2.mets)

Base.copy(model::LinearModel) =
    LinearModel(model.S, model.b, model.c, model.xl, model.xu, model.rxns, model.mets)

Base.isequal(model1::CoupledLinearModel, model2::CoupledLinearModel) =
    isequal(model1.lm, model2.lm) &&
    isequal(model1.C, model2.C) &&
    isequal(model1.cl, model2.cl) &&
    isequal(model1.cu, model2.cu)

Base.copy(model::CoupledLinearModel) =
    CoupledLinearModel(model.lm, model.C, model.cl, model.cu)

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

Return two dictionaries of metabolite `id`s mapped to reactions that consume or 
produce them given the flux distribution supplied in `fluxdict`.
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
The JuMP `set_normalized_rhs` function is a little confusing, 
so this function simplifies setting constraints. In short, JuMP
uses a normalized right hand side representation of constraints, 
which means that lower bounds have their sign flipped. This function
does this for you, so you don't have to remember to do this whenever you
change the constraints. 

Just supply the constraint `index` and the JuMP model (`opt_model`) that 
will be solved, and the variable's bounds will be set to `ub` and `lb`.
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
    change_constraint(reaction::Reaction, lb, ub)

Change the lower and upper bounds (`lb` and `ub` respectively) of `reaction`.
"""
function change_constraint(reaction::Reaction, lb, ub)
    (model, opt_model) -> begin
        ind = model.reactions[reaction]
        set_bound(ind, opt_model, lb = lb, ub = ub)
    end
end

"""
    change_solver_attribute(option_key, option_val)

Change a solver attribute. These attributes are solver specific,
refer the either JuMP or the solver you are using's documentation.
"""
function change_solver_attribute(option_key, option_val)
    (model, opt_model) -> begin
        COBREXA.JuMP.set_optimizer_attribute(opt_model, option_key, option_val)
    end
end

"""
    change_sense(objective_sense)

Change the objective sense. 
Possible arguments are `MOI.MAX_SENSE` and `MOI.MIN_SENSE`.
Note, [`change_objective`](@ref) sets the sense of the objective,
so it doesn't make sense to use this function AND [`change_objective`](@ref) simultaneously. 
"""
function change_sense(objective_sense)
    (model, opt_model) -> begin
        COBREXA.JuMP.set_objective_sense(opt_model, objective_sense)
    end
end

"""
    change_objective(objective_functions::Union{Reaction, Array{Reaction, 1}}; weights=[], sense=MOI.MAX_SENSE)

Callback function to change the objective function used in a constraint based analysis function. 
`objective_functions` can be a single reaction or an array of reactions (of type `Reaction`).
Optionally specify their `weights` and the sense of the new objective (`MOI.MAX_SENSE`, `MOI.MIN_SENSE`).
Note, this function sets the sense of the objective to `MOI.MAX_SENSE` by default if not specified.
"""
function change_objective(
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
    change_solver(optimizer)

Change the solver (`optimizer`) used to solve the model.
Typically the solver is specified as a required argument in a function. 
However, this function is useful if the problem has multiple subparts that 
require different solvers.

See also: [`parsimonious_flux_balance_analysis`](@ref)
"""
function change_solver(optimizer)
    (model, opt_model) -> begin
        COBREXA.JuMP.set_optimizer(opt_model, optimizer)
    end
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
    test_samples(samples::Array{Float64, 2}, mass_balance, balance, lbs, ubs)

Test if `samples` generated by some sampling algorithm, e.g. `hit_and_run`, are
consitent with the model's constraints. That is, if the `mass balances` and flux bound
constraints (`ubs`, `lbs`) are satisfied. Returns a vector of indices where the
constraints are violated. 

See also [`hit_and_run`](@ref), [`artificially_centered_hit_and_run`](@ref)
"""
function test_samples(samples::Array{Float64,2}, mass_balance, balance, lbs, ubs)
    violations = Int64[]
    tol = 1e-6
    for i = 1:size(samples, 2)
        if isapprox(mass_balance * samples[:, i], balance; atol = tol)
            equality = true
        else
            equality = false
        end
        inequality =
            all(abs.(lbs .- samples[:, i]) .<= tol) .==
            all(abs.(samples[:, i] .- ubs) .<= tol)
        if !all([equality, inequality])
            push!(violations, i)
        end
    end

    return violations
end
